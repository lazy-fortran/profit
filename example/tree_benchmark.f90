program tree_benchmark
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use omp_lib, only: omp_set_num_threads, omp_get_wtime
    use profit_status, only: profit_success
    use profit_models_tree, only: tree_options_t, tree_bins_t, regression_tree_t
    implicit none
    integer, parameter :: p = 16, n = 200000, repetitions = 5
    type(tree_options_t) :: options
    type(tree_bins_t) :: bins
    type(regression_tree_t) :: tree
    real(dp), allocatable :: x(:, :), y(:), w(:), pred(:), reference(:)
    integer, allocatable :: binned(:, :)
    logical, allocatable :: fit_mask(:)
    real(dp) :: started, elapsed_fit, elapsed_predict
    integer :: i, f, status, threads, rep

    allocate(x(p, n), y(n), w(n), pred(n), reference(n), binned(p, n), fit_mask(n))
    do i = 1, n
        do f = 1, p
            x(f, i) = sin(real(i*f, dp)/31) + cos(real(i + f, dp)/17)
        end do
        y(i) = x(2, i)**2 + x(3, i) - 2*x(4, i)*x(5, i)
    end do
    w = 1
    fit_mask = .true.
    options = tree_options_t(max_depth=3, min_leaf=100, max_bins=32, leaf_l2=1.0_dp)
    call bins%fit(x, fit_mask, options, status)
    if (status /= profit_success) error stop "bin fit failed"
    call bins%transform(x, binned, status)
    if (status /= profit_success) error stop "bin transform failed"
    write(*, "(a)") "threads,rows,features,repetitions,fit_seconds,predict_seconds,max_difference"
    do threads = 1, 4, 3
        call omp_set_num_threads(threads)
        started = omp_get_wtime()
        do rep = 1, repetitions
            call tree%fit(bins, binned, y, w, fit_mask, options, status)
            if (status /= profit_success) error stop "tree fit failed"
        end do
        elapsed_fit = (omp_get_wtime() - started)/repetitions
        started = omp_get_wtime()
        do rep = 1, repetitions
            call tree%predict(x, pred, status)
            if (status /= profit_success) error stop "tree prediction failed"
        end do
        elapsed_predict = (omp_get_wtime() - started)/repetitions
        if (threads == 1) reference = pred
        write(*, "(i0,3(',',i0),3(',',es20.12))") threads, n, p, repetitions, &
            elapsed_fit, elapsed_predict, maxval(abs(pred - reference))
    end do
end program tree_benchmark
