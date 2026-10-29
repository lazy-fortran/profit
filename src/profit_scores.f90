module profit_scores
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use profit_status
    use profit_generated_score_leaf, only: score_leaf
    use profit_generated_square_leaf, only: square_leaf
    use profit_generated_rmse_leaf, only: rmse_leaf
    implicit none
    private
    public :: gaussian_scores
contains
    subroutine gaussian_scores(y, mean, stddev, z, rmse, nll, coverage, status)
        real(dp), intent(in) :: y(:), mean(:), stddev(:), z
        real(dp), intent(out) :: rmse, nll, coverage
        integer, intent(out) :: status
        real(dp) :: residual, standardized, point_nll, scale, squared, loss, square
        integer :: i, n, covered, invalid
        rmse = 0.0_dp
        nll = 0.0_dp
        coverage = 0.0_dp
        status = profit_invalid_shape
        n = size(y)
        if (n == 0) return
        if (size(mean) /= n .or. size(stddev) /= n) return
        status = profit_invalid_input
        if (.not. ieee_is_finite(z)) return
        if (z < 0.0_dp) return
        if (.not. all(ieee_is_finite(y))) return
        if (.not. all(ieee_is_finite(mean))) return
        if (.not. all(ieee_is_finite(stddev))) return
        if (any(stddev <= 0.0_dp)) return
        scale = maxval(abs(y - mean))
        status = profit_range_error
        if (.not. ieee_is_finite(scale)) return
        squared = 0.0_dp
        loss = 0.0_dp
        covered = 0
        invalid = 0
        !$omp parallel do if(n > 4096) default(none) &
        !$omp shared(n,y,mean,stddev,z,scale) &
        !$omp private(i,residual,standardized,point_nll) &
        !$omp private(square) &
        !$omp reduction(+:squared,loss,covered,invalid)
        do i = 1, n
            residual = y(i) - mean(i)
            call score_leaf(residual, stddev(i), standardized, point_nll)
            if (.not. ieee_is_finite(point_nll)) invalid = invalid + 1
            if (scale > 0.0_dp) then
                call square_leaf(residual, scale, square)
                squared = squared + square
            end if
            loss = loss + point_nll/real(n, dp)
            if (abs(standardized) <= z) covered = covered + 1
        end do
        !$omp end parallel do
        if (invalid /= 0) return
        if (.not. ieee_is_finite(loss)) return
        call rmse_leaf(scale, squared, real(n, dp), rmse)
        nll = loss
        coverage = real(covered, dp)/real(n, dp)
        status = profit_success
    end subroutine gaussian_scores
end module profit_scores
