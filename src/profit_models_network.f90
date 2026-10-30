module profit_models_network
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use profit_status, only: profit_success, profit_invalid_shape, &
        profit_invalid_input, profit_range_error
    use profit_generated_network_activation, only: network_activation
    use profit_generated_network_positive, only: network_positive
    use profit_generated_network_small_log, only: network_small_log
    implicit none
    private
    public :: network_parameter_count, network_predict, network_jvp, network_vjp
    public :: network_positive_map
contains
    pure integer function network_parameter_count(ninput, nhidden) result(count)
        integer, intent(in) :: ninput, nhidden
        integer(int64) :: required
        count = 0
        if (ninput < 1 .or. nhidden < 0) return
        required = int(ninput, int64) + 1_int64
        if (nhidden > 0) required = int(nhidden, int64)* &
            (int(ninput, int64) + 2_int64) + 1_int64
        if (required > int(huge(count), int64)) return
        count = int(required)
    end function network_parameter_count

    subroutine network_predict(ninput, nhidden, theta, features, values, status)
        integer, intent(in) :: ninput, nhidden
        real(dp), intent(in) :: theta(:), features(:, :)
        real(dp), intent(out) :: values(:)
        integer, intent(out) :: status
        real(dp) :: z, activation, derivative, value
        integer :: i, j, first, bias, output, invalid, n
        values = 0.0_dp
        call validate(ninput, nhidden, theta, features, size(values), status)
        if (status /= profit_success) return
        n = size(values)
        bias = ninput*nhidden
        output = bias + nhidden
        invalid = 0
        !$omp parallel do if(n > 4096) default(none) &
        !$omp shared(n,ninput,nhidden,theta,features,values,bias,output) &
        !$omp private(i,j,first,z,activation,derivative,value) reduction(+:invalid)
        do i = 1, n
            if (nhidden == 0) then
                value = theta(ninput + 1) + dot_product(theta(1:ninput), features(:, i))
            else
                value = theta(size(theta))
                do j = 1, nhidden
                    first = (j - 1)*ninput
                    z = theta(bias + j) + &
                        dot_product(theta(first + 1:first + ninput), features(:, i))
                    if (.not. ieee_is_finite(z)) exit
                    call network_activation(z, activation, derivative)
                    value = value + theta(output + j)*activation
                end do
                if (.not. ieee_is_finite(z)) then
                    invalid = invalid + 1
                    cycle
                end if
            end if
            if (.not. ieee_is_finite(value)) then
                invalid = invalid + 1
                cycle
            end if
            values(i) = value
        end do
        !$omp end parallel do
        status = profit_range_error
        if (invalid /= 0) then
            values = 0.0_dp
            return
        end if
        status = profit_success
    end subroutine network_predict

    subroutine network_jvp(ninput, nhidden, theta, features, theta_direction, &
            feature_direction, dvalues, status)
        integer, intent(in) :: ninput, nhidden
        real(dp), intent(in) :: theta(:), features(:, :), theta_direction(:)
        real(dp), intent(in) :: feature_direction(:, :)
        real(dp), intent(out) :: dvalues(:)
        integer, intent(out) :: status
        real(dp) :: z, dz, activation, derivative, value
        integer :: i, j, first, bias, output, invalid, n
        dvalues = 0.0_dp
        call validate(ninput, nhidden, theta, features, size(dvalues), status)
        if (status /= profit_success) return
        status = profit_invalid_shape
        if (size(theta_direction) /= size(theta)) return
        if (any(shape(feature_direction) /= shape(features))) return
        status = profit_invalid_input
        if (.not. all(ieee_is_finite(theta_direction))) return
        if (.not. all(ieee_is_finite(feature_direction))) return
        n = size(dvalues)
        bias = ninput*nhidden
        output = bias + nhidden
        invalid = 0
        !$omp parallel do if(n > 4096) default(none) &
        !$omp shared(n,ninput,nhidden,theta,features,theta_direction) &
        !$omp shared(feature_direction,dvalues,bias,output) &
        !$omp private(i,j,first,z,dz,activation,derivative,value) reduction(+:invalid)
        do i = 1, n
            if (nhidden == 0) then
                value = theta_direction(ninput + 1) + &
                    dot_product(theta_direction(1:ninput), features(:, i)) + &
                    dot_product(theta(1:ninput), feature_direction(:, i))
            else
                value = theta_direction(size(theta))
                do j = 1, nhidden
                    first = (j - 1)*ninput
                    z = theta(bias + j) + &
                        dot_product(theta(first + 1:first + ninput), features(:, i))
                    dz = theta_direction(bias + j) + &
                        dot_product(theta_direction(first + 1:first + ninput), &
                        features(:, i)) + &
                        dot_product(theta(first + 1:first + ninput), &
                        feature_direction(:, i))
                    if (.not. ieee_is_finite(z)) exit
                    if (.not. ieee_is_finite(dz)) exit
                    call network_activation(z, activation, derivative)
                    value = value + theta_direction(output + j)*activation + &
                        theta(output + j)*derivative*dz
                end do
                if (.not. ieee_is_finite(z) .or. .not. ieee_is_finite(dz)) then
                    invalid = invalid + 1
                    cycle
                end if
            end if
            if (.not. ieee_is_finite(value)) then
                invalid = invalid + 1
                cycle
            end if
            dvalues(i) = value
        end do
        !$omp end parallel do
        status = profit_range_error
        if (invalid /= 0) then
            dvalues = 0.0_dp
            return
        end if
        status = profit_success
    end subroutine network_jvp

    subroutine network_vjp(ninput, nhidden, theta, features, output_adjoint, &
            theta_adjoint, feature_adjoint, status)
        integer, intent(in) :: ninput, nhidden
        real(dp), intent(in) :: theta(:), features(:, :), output_adjoint(:)
        real(dp), intent(out) :: theta_adjoint(:), feature_adjoint(:, :)
        integer, intent(out) :: status
        real(dp) :: z, activation, derivative, hidden_adjoint, weight
        integer :: i, j, k, first, bias, output, invalid, n
        theta_adjoint = 0.0_dp
        feature_adjoint = 0.0_dp
        call validate(ninput, nhidden, theta, features, size(output_adjoint), status)
        if (status /= profit_success) return
        status = profit_invalid_shape
        if (size(theta_adjoint) /= size(theta)) return
        if (any(shape(feature_adjoint) /= shape(features))) return
        status = profit_invalid_input
        if (.not. all(ieee_is_finite(output_adjoint))) return
        n = size(output_adjoint)
        bias = ninput*nhidden
        output = bias + nhidden
        invalid = 0
        !$omp parallel do if(n > 4096) default(none) &
        !$omp shared(n,ninput,nhidden,theta,features,output_adjoint,feature_adjoint) &
        !$omp shared(bias,output) private(i,j,k,first,z,activation,derivative) &
        !$omp private(hidden_adjoint,weight) reduction(+:invalid,theta_adjoint)
        do i = 1, n
            weight = output_adjoint(i)
            if (nhidden == 0) then
                theta_adjoint(ninput + 1) = theta_adjoint(ninput + 1) + weight
                do k = 1, ninput
                    theta_adjoint(k) = theta_adjoint(k) + weight*features(k, i)
                    feature_adjoint(k, i) = weight*theta(k)
                end do
            else
                theta_adjoint(size(theta)) = theta_adjoint(size(theta)) + weight
                do j = 1, nhidden
                    first = (j - 1)*ninput
                    z = theta(bias + j) + &
                        dot_product(theta(first + 1:first + ninput), features(:, i))
                    if (.not. ieee_is_finite(z)) exit
                    call network_activation(z, activation, derivative)
                    theta_adjoint(output + j) = theta_adjoint(output + j) + &
                        weight*activation
                    hidden_adjoint = weight*theta(output + j)*derivative
                    theta_adjoint(bias + j) = theta_adjoint(bias + j) + hidden_adjoint
                    do k = 1, ninput
                        theta_adjoint(first + k) = theta_adjoint(first + k) + &
                            hidden_adjoint*features(k, i)
                        feature_adjoint(k, i) = feature_adjoint(k, i) + &
                            hidden_adjoint*theta(first + k)
                    end do
                end do
                if (.not. ieee_is_finite(z)) invalid = invalid + 1
            end if
        end do
        !$omp end parallel do
        status = profit_range_error
        if (invalid /= 0 .or. .not. all(ieee_is_finite(theta_adjoint)) .or. &
            .not. all(ieee_is_finite(feature_adjoint))) then
            theta_adjoint = 0.0_dp
            feature_adjoint = 0.0_dp
            return
        end if
        status = profit_success
    end subroutine network_vjp

    subroutine network_positive_map(scores, values, derivatives, status)
        real(dp), intent(in) :: scores(:)
        real(dp), intent(out) :: values(:), derivatives(:)
        integer, intent(out) :: status
        real(dp) :: tail, correction, slope
        integer :: i, n
        values = 0.0_dp
        derivatives = 0.0_dp
        n = size(scores)
        status = profit_invalid_shape
        if (n == 0) return
        if (size(values) /= n .or. size(derivatives) /= n) return
        status = profit_invalid_input
        if (.not. all(ieee_is_finite(scores))) return
        !$omp parallel do if(n > 4096) default(none) &
        !$omp shared(n,scores,values,derivatives) private(i,tail,correction,slope)
        do i = 1, n
            call network_positive(abs(scores(i)), tail, correction, slope)
            if (tail < 1e-3_dp) call network_small_log(tail, correction)
            values(i) = correction
            if (scores(i) >= 0.0_dp) then
                values(i) = scores(i) + correction
                derivatives(i) = slope
            else
                derivatives(i) = tail*slope
            end if
        end do
        !$omp end parallel do
        status = profit_success
    end subroutine network_positive_map

    subroutine validate(ninput, nhidden, theta, features, n, status)
        integer, intent(in) :: ninput, nhidden, n
        real(dp), intent(in) :: theta(:), features(:, :)
        integer, intent(out) :: status
        integer :: required
        status = profit_invalid_shape
        required = network_parameter_count(ninput, nhidden)
        if (required == 0 .or. n == 0) return
        if (size(theta) /= required) return
        if (size(features, 1) /= ninput .or. size(features, 2) /= n) return
        status = profit_invalid_input
        if (.not. all(ieee_is_finite(theta))) return
        if (.not. all(ieee_is_finite(features))) return
        status = profit_success
    end subroutine validate
end module profit_models_network
