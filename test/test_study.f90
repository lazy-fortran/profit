module study_tests
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use profit_scores, only: gaussian_scores
    use profit_replicates, only: replicate_summary
    use profit_design, only: expected_improvement
    use profit_output, only: write_study_csv
    use profit_status
    implicit none
contains
    subroutine collect(tests)
        type(unittest_type), allocatable, intent(out) :: tests(:)
        allocate(tests(7))
        tests(1) = new_unittest("Gaussian scores and units", test_scores)
        tests(2) = new_unittest("independent replicate uncertainty", test_replicates)
        tests(3) = new_unittest("EI exact cases and quadrature", test_ei)
        tests(4) = new_unittest("EI derivative products", test_derivatives)
        tests(5) = new_unittest("invalid inputs", test_invalid)
        tests(6) = new_unittest("parallel batched operations", test_batches)
        tests(7) = new_unittest("CSV records and provenance", test_csv)
    end subroutine collect

    subroutine test_scores(error)
        type(error_type), allocatable, intent(out) :: error
        real(dp) :: rmse, nll, coverage, scaled_rmse, scaled_nll, scaled_coverage
        integer :: status
        call gaussian_scores([0.0_dp, 2.0_dp], [0.0_dp, 0.0_dp], &
            [1.0_dp, 1.0_dp], 1.0_dp, rmse, nll, coverage, status)
        call check(error, status, profit_success)
        if (allocated(error)) return
        call check(error, rmse, sqrt(2.0_dp), thr=1e-14_dp)
        if (allocated(error)) return
        call check(error, nll, log(2*acos(-1.0_dp))/2 + 1, thr=1e-14_dp)
        if (allocated(error)) return
        call check(error, coverage, 0.5_dp)
        if (allocated(error)) return
        call gaussian_scores([0.0_dp, 200.0_dp], [0.0_dp, 0.0_dp], &
            [100.0_dp, 100.0_dp], 1.0_dp, scaled_rmse, scaled_nll, &
            scaled_coverage, status)
        call check(error, status, profit_success)
        if (allocated(error)) return
        call check(error, scaled_rmse, 100*rmse, thr=1e-12_dp)
        if (allocated(error)) return
        call check(error, scaled_nll, nll + log(100.0_dp), thr=1e-14_dp)
        if (allocated(error)) return
        call check(error, scaled_coverage, coverage)
    end subroutine test_scores

    subroutine test_replicates(error)
        type(error_type), allocatable, intent(out) :: error
        real(dp) :: mean, se, lower, upper
        integer :: status
        call replicate_summary([1.0_dp, 2.0_dp, 3.0_dp], 2.0_dp, &
            mean, se, lower, upper, status)
        call check(error, status, profit_success)
        if (allocated(error)) return
        call check(error, mean, 2.0_dp, thr=1e-14_dp)
        if (allocated(error)) return
        call check(error, se, 1/sqrt(3.0_dp), thr=1e-14_dp)
        if (allocated(error)) return
        call check(error, lower, 2 - 2/sqrt(3.0_dp), thr=1e-14_dp)
        if (allocated(error)) return
        call check(error, upper, 2 + 2/sqrt(3.0_dp), thr=1e-14_dp)
        if (allocated(error)) return
        call replicate_summary([4.0_dp, 4.0_dp], 2.0_dp, &
            mean, se, lower, upper, status)
        call check(error, status == profit_success .and. se == 0.0_dp)
    end subroutine test_replicates

    subroutine test_ei(error)
        type(error_type), allocatable, intent(out) :: error
        real(dp), parameter :: mus(5) = [-4.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 4.0_dp]
        real(dp) :: ei(5), reflected(5), exact, z, weight, integral, step, tail
        integer :: status, i, j
        call expected_improvement(mus, [1.0_dp, 1.0_dp, 1.0_dp, &
            1.0_dp, 1.0_dp], 0.0_dp, .false., ei, status)
        call check(error, status, profit_success)
        if (allocated(error)) return
        call check(error, ei(3), 1/sqrt(2*acos(-1.0_dp)), thr=1e-14_dp)
        if (allocated(error)) return
        call expected_improvement(-mus, [1.0_dp, 1.0_dp, 1.0_dp, &
            1.0_dp, 1.0_dp], 0.0_dp, .true., reflected, status)
        call check(error, maxval(abs(ei - reflected)) < 1e-14_dp)
        if (allocated(error)) return
        ! Integrate positive improvement against a standard normal density.
        do i = 1, size(mus)
            tail = -mus(i)
            step = (12.0_dp - tail)/20000
            integral = 0.0_dp
            do j = 0, 20000
                z = tail + j*step
                weight = 2.0_dp
                if (mod(j, 2) == 1) weight = 4.0_dp
                if (j == 0 .or. j == 20000) weight = 1.0_dp
                integral = integral + weight*max(mus(i) + z, 0.0_dp)* &
                    exp(-z*z/2)/sqrt(2*acos(-1.0_dp))
            end do
            exact = integral*step/3
            call check(error, ei(i), exact, thr=2e-12_dp)
            if (allocated(error)) return
        end do
        call expected_improvement([2.0_dp, -2.0_dp, 0.0_dp, &
            -100.0_dp, 100.0_dp], [0.0_dp, 0.0_dp, 0.0_dp, &
            1.0_dp, 1.0_dp], 0.0_dp, .false., ei, status)
        call check(error, status, profit_success)
        if (allocated(error)) return
        call check(error, all(ei == [2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 100.0_dp]))
    end subroutine test_ei

    subroutine test_derivatives(error)
        type(error_type), allocatable, intent(out) :: error
        real(dp), parameter :: h = 1e-5_dp
        real(dp) :: ei(1), dm(1), ds(1), plus(1), minus(1)
        integer :: status, direction
        do direction = 0, 1
            call expected_improvement([0.3_dp], [0.7_dp], 0.2_dp, &
                direction == 0, ei, status, dm, ds)
            call check(error, status, profit_success)
            if (allocated(error)) return
            call expected_improvement([0.3_dp + h], [0.7_dp], 0.2_dp, &
                direction == 0, plus, status)
            call expected_improvement([0.3_dp - h], [0.7_dp], 0.2_dp, &
                direction == 0, minus, status)
            call check(error, dm(1), (plus(1) - minus(1))/(2*h), thr=2e-10_dp)
            if (allocated(error)) return
            call expected_improvement([0.3_dp], [0.7_dp + h], 0.2_dp, &
                direction == 0, plus, status)
            call expected_improvement([0.3_dp], [0.7_dp - h], 0.2_dp, &
                direction == 0, minus, status)
            call check(error, ds(1), (plus(1) - minus(1))/(2*h), thr=2e-10_dp)
            if (allocated(error)) return
        end do
    end subroutine test_derivatives

    subroutine test_invalid(error)
        type(error_type), allocatable, intent(out) :: error
        real(dp) :: rmse, nll, coverage, mean, se, lower, upper, ei(1), nan_input(1)
        integer :: status
        call gaussian_scores([0.0_dp], [0.0_dp], [0.0_dp], 1.0_dp, &
            rmse, nll, coverage, status)
        call check(error, status, profit_invalid_input)
        if (allocated(error)) return
        call gaussian_scores([0.0_dp], [0.0_dp, 1.0_dp], [1.0_dp], 1.0_dp, &
            rmse, nll, coverage, status)
        call check(error, status, profit_invalid_shape)
        if (allocated(error)) return
        call replicate_summary([1.0_dp], 2.0_dp, mean, se, lower, upper, status)
        call check(error, status, profit_invalid_shape)
        if (allocated(error)) return
        nan_input = ieee_value(0.0_dp, ieee_quiet_nan)
        call expected_improvement(nan_input, [1.0_dp], 0.0_dp, .true., ei, status)
        call check(error, status, profit_invalid_input)
        if (allocated(error)) return
        call expected_improvement([0.0_dp], [-1.0_dp], 0.0_dp, .true., ei, status)
        call check(error, status, profit_invalid_input)
    end subroutine test_invalid

    subroutine test_batches(error)
        type(error_type), allocatable, intent(out) :: error
        integer, parameter :: n = 20000
        real(dp) :: y(n), mu(n), sd(n), ei(n), dm(n), ds(n)
        real(dp) :: rmse, nll, coverage, average, se, lower, upper
        integer :: status
        y(1:n:2) = 1.0_dp
        y(2:n:2) = 3.0_dp
        mu = 2.0_dp
        sd = 1.0_dp
        call gaussian_scores(y, mu, sd, 1.0_dp, rmse, nll, coverage, status)
        call check(error, status, profit_success)
        if (allocated(error)) return
        call check(error, rmse, 1.0_dp, thr=1e-13_dp)
        if (allocated(error)) return
        call check(error, nll, (log(2*acos(-1.0_dp)) + 1)/2, thr=1e-12_dp)
        if (allocated(error)) return
        call check(error, coverage, 1.0_dp)
        if (allocated(error)) return
        call replicate_summary(y, 2.0_dp, average, se, lower, upper, status)
        call check(error, status, profit_success)
        if (allocated(error)) return
        call check(error, se, 1/sqrt(real(n - 1, dp)), thr=1e-14_dp)
        if (allocated(error)) return
        call expected_improvement(mu, sd, 2.0_dp, .true., ei, status, dm, ds)
        call check(error, status, profit_success)
        if (allocated(error)) return
        call check(error, maxval(abs(ei - 1/sqrt(2*acos(-1.0_dp)))) < 1e-14_dp)
        if (allocated(error)) return
        call check(error, maxval(abs(dm + 0.5_dp)) < 1e-14_dp)
    end subroutine test_batches

    subroutine test_csv(error)
        type(error_type), allocatable, intent(out) :: error
        character(*), parameter :: path = "profit-test-study.csv"
        character(512) :: line
        character(6), parameter :: keys(1) = ["source"]
        character(17), parameter :: metadata(1) = ['synthetic,"known"']
        character(8), parameter :: columns(1) = ["rmse [K]"]
        real(dp), parameter :: table(1, 1) = reshape([1.25_dp], [1, 1])
        real(dp) :: number
        integer :: status, unit, ios
        call write_study_csv(path, keys, metadata, columns, table, status)
        call check(error, status, profit_success)
        if (allocated(error)) return
        open (newunit=unit, file=path, status="old", iostat=ios)
        call check(error, ios, 0)
        if (allocated(error)) return
        read (unit, '(a)', iostat=ios) line
        call check(error, trim(line), &
            'record_type,metadata_key,metadata_value,"rmse [K]"')
        if (.not. allocated(error)) then
            read (unit, '(a)', iostat=ios) line
            call check(error, trim(line), 'metadata,"source","synthetic,""known""",')
        end if
        if (.not. allocated(error)) then
            read (unit, '(a)', iostat=ios) line
            read (line(8:), *, iostat=ios) number
            call check(error, number, 1.25_dp)
        end if
        close (unit, status="delete")
    end subroutine test_csv
end module study_tests

program test_study
    use, intrinsic :: iso_fortran_env, only: error_unit
    use testdrive, only: run_testsuite
    use study_tests, only: collect
    implicit none
    integer :: failures
    failures = 0
    call run_testsuite(collect, error_unit, failures)
    if (failures /= 0) error stop 1
end program test_study
