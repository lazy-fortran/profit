program expected_improvement_curve
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use profit_design, only: expected_improvement
    use profit_status, only: profit_success
    implicit none
    integer, parameter :: points = 121
    real(dp) :: mean(points), stddev(points), analytic(points), reference(points)
    real(dp) :: dmean(points), dstddev(points), error(points)
    character(4096) :: path
    integer :: i, unit, status, ios

    path = "expected_improvement_curve.csv"
    if (command_argument_count() > 1) error stop "supply at most one CSV path"
    if (command_argument_count() == 1) call get_command_argument(1, path)
    do i = 1, points
        mean(i) = -3.0_dp + 6.0_dp*real(i - 1, dp)/real(points - 1, dp)
    end do
    stddev = 1.0_dp
    call expected_improvement(mean, stddev, 0.0_dp, .false., analytic, &
        status, dmean, dstddev)
    if (status /= profit_success) error stop "expected improvement failed"

    ! A numerical reference integrates positive improvement against the density.
    !$omp parallel do default(none) shared(mean,reference) private(i)
    do i = 1, points
        reference(i) = quadrature_reference(mean(i))
    end do
    !$omp end parallel do
    error = abs(analytic - reference)
    if (maxval(error) > 1e-11_dp) error stop "quadrature comparison failed"

    open(newunit=unit, file=trim(path), status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "cannot open curve CSV"
    write(unit, '(a)', iostat=ios) &
        "mean_minus_incumbent,analytic_ei,quadrature_ei,absolute_error,dmean,dstddev"
    if (ios /= 0) error stop "cannot write CSV header"
    do i = 1, points
        write(unit, '(*(es25.17e3,:,","))', iostat=ios) mean(i), analytic(i), &
            reference(i), error(i), dmean(i), dstddev(i)
        if (ios /= 0) error stop "cannot write curve CSV"
    end do
    close(unit, iostat=ios)
    if (ios /= 0) error stop "cannot close curve CSV"
    print '(a,es12.4)', "Maximum absolute quadrature difference: ", maxval(error)

contains

    pure function quadrature_reference(mu) result(value)
        real(dp), intent(in) :: mu
        real(dp) :: value, step, z, weight, integral
        integer, parameter :: intervals = 20000
        integer :: j
        ! At sigma=1, improvement is mu+z; its positive region starts at -mu.
        step = (12.0_dp + mu)/real(intervals, dp)
        integral = 0.0_dp
        do j = 0, intervals
            z = -mu + real(j, dp)*step
            weight = 2.0_dp
            if (mod(j, 2) == 1) weight = 4.0_dp
            if (j == 0 .or. j == intervals) weight = 1.0_dp
            integral = integral + weight*max(mu + z, 0.0_dp)* &
                exp(-z*z/2)/sqrt(2*acos(-1.0_dp))
        end do
        value = integral*step/3
    end function quadrature_reference
end program expected_improvement_curve
