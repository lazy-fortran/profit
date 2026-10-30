module tree_tests
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use omp_lib, only: omp_set_num_threads, omp_get_max_threads
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use profit_status, only: profit_success, profit_invalid_input, &
        profit_invalid_shape, profit_io_error
    use profit_models_tree, only: tree_options_t, tree_bins_t, regression_tree_t
    implicit none
contains
    subroutine collect(tests)
        type(unittest_type), allocatable, intent(out) :: tests(:)
        allocate(tests(7))
        tests(1) = new_unittest("weighted CART versus exhaustive split", test_oracle)
        tests(2) = new_unittest("FIT-only bins and targets", test_fit_only)
        tests(3) = new_unittest("constant response and ridge leaves", test_constant)
        tests(4) = new_unittest("leaf editing and serialization", test_leaves)
        tests(5) = new_unittest("invalid tree inputs", test_invalid)
        tests(6) = new_unittest("parallel histogram determinism", test_parallel)
        tests(7) = new_unittest("recurrent objective descent", test_recurrent)
    end subroutine collect

    subroutine test_oracle(error)
        type(error_type), allocatable, intent(out) :: error
        type(tree_options_t) :: options
        type(tree_bins_t) :: bins
        type(regression_tree_t) :: tree
        real(dp) :: x(2, 12), y(12), w(12), pred(12), best(12), candidate(12)
        real(dp) :: loss, best_loss, cutoff, a, b
        logical :: mask(12), left(12)
        integer :: i, f, j, status, binned(2, 12)

        do i = 1, 12
            x(1, i) = real(mod(5*i, 13), dp)
            x(2, i) = real(mod(7*i, 13), dp)
            y(i) = 0.2_dp*i
            if (x(2, i) > 5) y(i) = y(i) + 4
            w(i) = 1 + real(mod(i, 3), dp)
        end do
        mask = .true.
        options = tree_options_t(max_depth=1, min_leaf=2, max_bins=16, leaf_l2=0.3_dp)
        call bins%fit(x, mask, options, status)
        call check(error, status, profit_success)
        if (allocated(error)) return
        call bins%transform(x, binned, status)
        call check(error, status, profit_success)
        if (allocated(error)) return
        call tree%fit(bins, binned, y, w, mask, options, status)
        call check(error, status, profit_success)
        if (allocated(error)) return
        call tree%predict(x, pred, status)
        call check(error, status, profit_success)
        if (allocated(error)) return
        best_loss = huge(1.0_dp)
        do f = 1, 2
            do j = 1, 12
                cutoff = x(f, j)
                left = x(f, :) <= cutoff
                if (count(left) < 2 .or. count(.not. left) < 2) cycle
                a = sum(w*y, mask=left)/(sum(w, mask=left) + options%leaf_l2)
                b = sum(w*y, mask=.not. left)/(sum(w, mask=.not. left) + options%leaf_l2)
                candidate = merge(a, b, left)
                loss = sum(w*(candidate - y)**2) + options%leaf_l2*(a*a + b*b)
                if (loss < best_loss) then
                    best_loss = loss
                    best = candidate
                end if
            end do
        end do
        call check(error, maxval(abs(pred - best)) < 1e-13_dp)
        if (allocated(error)) return
        call check(error, tree%n_leaves(), 2)
    end subroutine test_oracle

    subroutine test_fit_only(error)
        type(error_type), allocatable, intent(out) :: error
        type(tree_options_t) :: options
        type(tree_bins_t) :: bins, other_bins
        type(regression_tree_t) :: tree, other_tree
        real(dp) :: x(1, 8), changed_x(1, 8), y(8), changed_y(8), w(8)
        real(dp) :: pred(6), other_pred(6), nan
        logical :: mask(8)
        integer :: binned(1, 8), other_binned(1, 8), status

        nan = ieee_value(0.0_dp, ieee_quiet_nan)
        x(1, :6) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]
        x(1, 7:) = nan
        y(:3) = 0
        y(4:6) = 3
        y(7:) = nan
        w = 1
        mask = [.true., .true., .true., .true., .true., .true., .false., .false.]
        changed_x = x
        changed_x(1, 7:) = [-10000.0_dp, 10000.0_dp]
        changed_y = y
        changed_y(7:) = [-1e9_dp, 1e9_dp]
        options = tree_options_t(max_depth=1, min_leaf=1, max_bins=16, leaf_l2=0.0_dp)
        call bins%fit(x, mask, options, status)
        call check(error, status, profit_success)
        if (allocated(error)) return
        call other_bins%fit(changed_x, mask, options, status)
        call check(error, status, profit_success)
        if (allocated(error)) return
        call bins%transform(x, binned, status, mask)
        call check(error, status, profit_success)
        if (allocated(error)) return
        call other_bins%transform(changed_x, other_binned, status, mask)
        call check(error, status, profit_success)
        if (allocated(error)) return
        call check(error, all(binned == other_binned))
        if (allocated(error)) return
        call tree%fit(bins, binned, y, w, mask, options, status)
        call check(error, status, profit_success)
        if (allocated(error)) return
        call other_tree%fit(other_bins, other_binned, changed_y, w, mask, options, status)
        call check(error, status, profit_success)
        if (allocated(error)) return
        call tree%predict(x(:, :6), pred, status)
        call other_tree%predict(x(:, :6), other_pred, status)
        call check(error, all(pred == other_pred))
        if (allocated(error)) return
        call check(error, all(pred == y(:6)))
    end subroutine test_fit_only

    subroutine test_constant(error)
        type(error_type), allocatable, intent(out) :: error
        type(tree_options_t) :: options
        type(tree_bins_t) :: bins
        type(regression_tree_t) :: tree
        real(dp) :: x(1, 4), y(4), w(4), pred(4)
        integer :: binned(1, 4), status
        logical :: mask(4)

        x = 1
        y = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
        w = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
        mask = .true.
        options = tree_options_t(max_depth=3, min_leaf=1, max_bins=8, leaf_l2=2.0_dp)
        call bins%fit(x, mask, options, status)
        call bins%transform(x, binned, status)
        call tree%fit(bins, binned, y, w, mask, options, status)
        call check(error, status, profit_success)
        if (allocated(error)) return
        call check(error, tree%n_leaves(), 1)
        if (allocated(error)) return
        call tree%predict(x, pred, status)
        call check(error, maxval(abs(pred - 2.5_dp)) < 1e-14_dp)
        if (allocated(error)) return
        x(1, :) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
        y = 0
        call bins%fit(x, mask, options, status)
        call bins%transform(x, binned, status)
        call tree%fit(bins, binned, y, w, mask, options, status)
        call check(error, status, profit_success)
        if (allocated(error)) return
        call check(error, tree%n_leaves(), 1)
    end subroutine test_constant

    subroutine test_leaves(error)
        type(error_type), allocatable, intent(out) :: error
        type(tree_options_t) :: options
        type(tree_bins_t) :: bins
        type(regression_tree_t) :: tree, loaded
        real(dp) :: x(1, 8), y(8), w(8), pred(8), reloaded(8)
        real(dp), allocatable :: leaves(:)
        logical :: mask(8)
        integer :: status, i, unit, binned(1, 8), indices(8)

        x(1, :) = [(real(i, dp), i = 1, 8)]
        y = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp]
        w = 1
        mask = .true.
        options = tree_options_t(max_depth=2, min_leaf=1, max_bins=16, leaf_l2=0.0_dp)
        call bins%fit(x, mask, options, status)
        call bins%transform(x, binned, status)
        call tree%fit(bins, binned, y, w, mask, options, status)
        call check(error, status, profit_success)
        if (allocated(error)) return
        call tree%leaf_indices(x, indices, status)
        call check(error, status, profit_success)
        if (allocated(error)) return
        allocate(leaves(tree%n_leaves()))
        call tree%get_leaf_values(leaves, status)
        call check(error, status, profit_success)
        if (allocated(error)) return
        do i = 1, size(leaves)
            leaves(i) = leaves(i) + real(i, dp)/7
        end do
        call tree%set_leaf_values(leaves, status)
        call check(error, status, profit_success)
        if (allocated(error)) return
        call tree%predict(x, pred, status)
        call check(error, all(pred == leaves(indices)))
        if (allocated(error)) return
        call tree%save("tree-roundtrip.bin", status)
        call check(error, status, profit_success)
        if (allocated(error)) return
        call loaded%load("tree-roundtrip.bin", status)
        call check(error, status, profit_success)
        if (allocated(error)) return
        call loaded%predict(x, reloaded, status)
        call check(error, all(pred == reloaded))
        if (allocated(error)) return
        open(newunit=unit, file="tree-roundtrip.bin", status="old")
        close(unit, status="delete")
        open(newunit=unit, file="tree-invalid.bin", status="replace")
        write(unit, "(a)") "invalid tree"
        close(unit)
        call loaded%load("tree-invalid.bin", status)
        call check(error, status, profit_io_error)
        open(newunit=unit, file="tree-invalid.bin", status="old")
        close(unit, status="delete")
    end subroutine test_leaves

    subroutine test_invalid(error)
        type(error_type), allocatable, intent(out) :: error
        type(tree_options_t) :: options
        type(tree_bins_t) :: bins
        type(regression_tree_t) :: tree
        real(dp) :: x(1, 2), y(2), w(2), pred(2), nan
        integer :: status, binned(1, 2), wrong_binned(2, 2)
        logical :: mask(2)

        nan = ieee_value(0.0_dp, ieee_quiet_nan)
        x(1, :) = [0.0_dp, 1.0_dp]
        y = 1
        w = 1
        mask = .true.
        options = tree_options_t(max_depth=1, min_leaf=1, max_bins=8, leaf_l2=0.0_dp)
        call bins%fit(x, [.true.], options, status)
        call check(error, status, profit_invalid_shape)
        if (allocated(error)) return
        options%max_depth = 11
        call bins%fit(x, mask, options, status)
        call check(error, status, profit_invalid_input)
        if (allocated(error)) return
        options%max_depth = 1
        call bins%fit(x, mask, options, status)
        call bins%transform(x, wrong_binned, status)
        call check(error, status, profit_invalid_shape)
        if (allocated(error)) return
        call bins%transform(x, binned, status)
        w(2) = -1
        call tree%fit(bins, binned, y, w, mask, options, status)
        call check(error, status, profit_invalid_input)
        if (allocated(error)) return
        w = 1
        y(2) = nan
        call tree%fit(bins, binned, y, w, mask, options, status)
        call check(error, status, profit_invalid_input)
        if (allocated(error)) return
        y = 1
        call tree%fit(bins, binned, y, w, mask, options, status)
        call check(error, status, profit_success)
        if (allocated(error)) return
        x(1, 2) = nan
        call tree%predict(x, pred, status)
        call check(error, status, profit_invalid_input)
        if (allocated(error)) return
        call check(error, all(pred == 0))
    end subroutine test_invalid

    subroutine test_parallel(error)
        type(error_type), allocatable, intent(out) :: error
        integer, parameter :: n = 4096, p = 5
        type(tree_options_t) :: options
        type(tree_bins_t) :: bins
        type(regression_tree_t) :: serial, parallel
        real(dp) :: x(p, n), y(n), w(n), pred(n), other(n)
        logical :: mask(n)
        integer :: i, f, status, previous_threads, binned(p, n)

        previous_threads = omp_get_max_threads()
        do i = 1, n
            do f = 1, p
                x(f, i) = sin(real(i*f, dp)/31)
            end do
            y(i) = x(2, i)**2 + x(3, i) - 2*x(4, i)*x(5, i)
            w(i) = 1 + real(mod(i, 7), dp)/7
        end do
        mask = .true.
        options = tree_options_t(max_depth=4, min_leaf=20, max_bins=32, leaf_l2=0.1_dp)
        call omp_set_num_threads(1)
        call bins%fit(x, mask, options, status)
        call bins%transform(x, binned, status)
        call serial%fit(bins, binned, y, w, mask, options, status)
        call check(error, status, profit_success)
        if (allocated(error)) return
        call serial%predict(x, pred, status)
        call omp_set_num_threads(4)
        call parallel%fit(bins, binned, y, w, mask, options, status)
        call check(error, status, profit_success)
        if (allocated(error)) return
        call parallel%predict(x, other, status)
        call omp_set_num_threads(previous_threads)
        call check(error, status, profit_success)
        if (allocated(error)) return
        call check(error, all(pred == other))
    end subroutine test_parallel

    subroutine test_recurrent(error)
        type(error_type), allocatable, intent(out) :: error
        integer, parameter :: n = 96
        type(tree_options_t) :: options
        type(tree_bins_t) :: bins
        type(regression_tree_t) :: tree
        real(dp) :: x(1, n), target(n), score(n), gradient(n), direction(n), w(n)
        real(dp) :: pseudo(n), trial_score(n)
        real(dp) :: state, loss, next_loss, step, derivative
        logical :: mask(n)
        integer :: i, trial, status, binned(1, n)

        state = 0
        do i = 1, n
            x(1, i) = sin(real(i, dp)/7)
            state = 0.8_dp*state + 0.2_dp*merge(2.0_dp, -1.0_dp, x(1, i) > 0)
            target(i) = state
        end do
        score = 0
        w = 1
        mask = .true.
        options = tree_options_t(max_depth=2, min_leaf=4, max_bins=16, leaf_l2=0.1_dp)
        call recurrence(score, target, loss, gradient)
        call bins%fit(x, mask, options, status)
        call bins%transform(x, binned, status)
        pseudo = -gradient
        call tree%fit(bins, binned, pseudo, w, mask, options, status)
        call check(error, status, profit_success)
        if (allocated(error)) return
        call tree%predict(x, direction, status)
        derivative = dot_product(gradient, direction)
        call check(error, derivative < 0)
        if (allocated(error)) return
        step = 1
        do trial = 1, 30
            trial_score = score + step*direction
            call recurrence(trial_score, target, next_loss, gradient)
            if (next_loss <= loss + 1e-4_dp*step*derivative) exit
            step = step/2
        end do
        call check(error, next_loss < loss)
        if (allocated(error)) return
        call check(error, trial <= 30)
    end subroutine test_recurrent

    subroutine recurrence(score, target, loss, gradient)
        real(dp), intent(in) :: score(:), target(:)
        real(dp), intent(out) :: loss, gradient(:)
        real(dp) :: residual(size(score)), state, adjoint
        integer :: i

        state = 0
        do i = 1, size(score)
            state = 0.8_dp*state + 0.2_dp*score(i)
            residual(i) = state - target(i)
        end do
        loss = sum(residual**2)/2
        adjoint = 0
        do i = size(score), 1, -1
            adjoint = residual(i) + 0.8_dp*adjoint
            gradient(i) = 0.2_dp*adjoint
        end do
    end subroutine recurrence
end module tree_tests

program test_tree
    use, intrinsic :: iso_fortran_env, only: error_unit
    use testdrive, only: run_testsuite
    use tree_tests, only: collect
    implicit none
    integer :: stat

    stat = 0
    call run_testsuite(collect, error_unit, stat)
    if (stat /= 0) error stop 1
end program test_tree
