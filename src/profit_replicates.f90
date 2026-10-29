module profit_replicates
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use profit_status
    use profit_generated_square_leaf, only: square_leaf
    use profit_generated_replicate_leaf, only: replicate_leaf
    implicit none
    private
    public :: replicate_summary
contains
    subroutine replicate_summary(values, critical, mean, se, lower, upper, status)
        real(dp), intent(in) :: values(:), critical
        real(dp), intent(out) :: mean, se, lower, upper
        integer, intent(out) :: status
        real(dp) :: average, scale, squared, deviation, uncertainty, left, right
        integer :: i, n
        mean = 0.0_dp
        se = 0.0_dp
        lower = 0.0_dp
        upper = 0.0_dp
        status = profit_invalid_shape
        n = size(values)
        if (n < 2) return
        status = profit_invalid_input
        if (.not. ieee_is_finite(critical)) return
        if (critical <= 0.0_dp) return
        if (.not. all(ieee_is_finite(values))) return
        scale = maxval(abs(values))
        average = 0.0_dp
        if (scale > 0.0_dp) average = scale*sum(values/scale/real(n, dp))
        status = profit_range_error
        if (.not. ieee_is_finite(average)) return
        scale = maxval(abs(values - average))
        if (.not. ieee_is_finite(scale)) return
        squared = 0.0_dp
        !$omp parallel do if(n > 4096) default(none) &
        !$omp shared(n,values,average,scale) private(i,deviation) reduction(+:squared)
        do i = 1, n
            if (scale == 0.0_dp) cycle
            call square_leaf(values(i) - average, scale, deviation)
            squared = squared + deviation
        end do
        !$omp end parallel do
        call replicate_leaf(average, scale, squared, real(n, dp), critical, &
            uncertainty, left, right)
        if (.not. ieee_is_finite(left) .or. .not. ieee_is_finite(right)) return
        mean = average
        se = uncertainty
        lower = left
        upper = right
        status = profit_success
    end subroutine replicate_summary
end module profit_replicates
