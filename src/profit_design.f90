module profit_design
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use profit_status
    use profit_generated_ei_leaf, only: ei_leaf
    implicit none
    private
    public :: expected_improvement
contains
    subroutine expected_improvement(mean, stddev, incumbent, minimize, ei, &
            status, dmean, dstddev)
        real(dp), intent(in) :: mean(:), stddev(:), incumbent
        logical, intent(in) :: minimize
        real(dp), intent(out) :: ei(:)
        integer, intent(out) :: status
        real(dp), intent(out), optional :: dmean(:), dstddev(:)
        real(dp) :: sign, improvement, probability, density, value
        integer :: i, n, invalid
        ei = 0.0_dp
        if (present(dmean)) dmean = 0.0_dp
        if (present(dstddev)) dstddev = 0.0_dp
        n = size(mean)
        status = profit_invalid_shape
        if (n == 0) return
        if (size(stddev) /= n .or. size(ei) /= n) return
        if (present(dmean)) then
            if (size(dmean) /= n) return
        end if
        if (present(dstddev)) then
            if (size(dstddev) /= n) return
        end if
        status = profit_invalid_input
        if (.not. ieee_is_finite(incumbent)) return
        if (.not. all(ieee_is_finite(mean))) return
        if (.not. all(ieee_is_finite(stddev))) return
        if (any(stddev < 0.0_dp)) return
        sign = 1.0_dp
        if (minimize) sign = -1.0_dp
        invalid = 0
        !$omp parallel do if(n > 4096) default(none) &
        !$omp shared(n,mean,stddev,incumbent,sign,ei,dmean,dstddev) &
        !$omp private(i,improvement,probability,density,value) reduction(+:invalid)
        do i = 1, n
            improvement = sign*(mean(i) - incumbent)
            if (.not. ieee_is_finite(improvement)) then
                invalid = invalid + 1
                cycle
            end if
            if (stddev(i) == 0.0_dp) then
                value = max(improvement, 0.0_dp)
                probability = 0.0_dp
                if (improvement > 0.0_dp) probability = 1.0_dp
                density = 0.0_dp
                if (improvement == 0.0_dp) then
                    call ei_leaf(0.0_dp, 1.0_dp, value, probability, density)
                    value = 0.0_dp
                    probability = 0.0_dp
                end if
            else if (abs(improvement)/40.0_dp > stddev(i)) then
                value = max(improvement, 0.0_dp)
                probability = 0.0_dp
                if (improvement > 0.0_dp) probability = 1.0_dp
                density = 0.0_dp
            else
                call ei_leaf(improvement, stddev(i), value, probability, density)
            end if
            if (.not. ieee_is_finite(value)) then
                invalid = invalid + 1
                cycle
            end if
            ei(i) = max(value, 0.0_dp)
            if (present(dmean)) dmean(i) = sign*probability
            if (present(dstddev)) dstddev(i) = density
        end do
        !$omp end parallel do
        status = profit_range_error
        if (invalid /= 0) then
            ei = 0.0_dp
            if (present(dmean)) dmean = 0.0_dp
            if (present(dstddev)) dstddev = 0.0_dp
            return
        end if
        status = profit_success
    end subroutine expected_improvement
end module profit_design
