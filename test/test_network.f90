module network_tests
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: &
        ieee_value, ieee_quiet_nan, ieee_is_finite
    use omp_lib, only: omp_set_num_threads
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use profit_models_network, only: network_parameter_count, network_predict, &
        network_jvp, network_vjp, network_positive_map
    use profit_status
    implicit none
contains
    subroutine collect(tests)
        type(unittest_type), allocatable, intent(out) :: tests(:)
        allocate (tests(6))
        tests(1) = new_unittest("affine and tanh values", test_values)
        tests(2) = new_unittest("parameter and input derivative products", &
            test_products)
        tests(3) = new_unittest("local covariance propagation", test_covariance)
        tests(4) = new_unittest("stable positive transformation", test_positive)
        tests(5) = new_unittest("invalid inputs and numerical range", test_invalid)
        tests(6) = new_unittest("parallel derivative agreement", test_parallel)
    end subroutine collect

    subroutine test_values(error)
        type(error_type), allocatable, intent(out) :: error
        real(dp) :: features(2, 3), values(3), theta(9), expected(3)
        integer :: i, status
        features = reshape([1.0_dp, -1.0_dp, 0.5_dp, 2.0_dp, -2.0_dp, &
            0.75_dp], shape(features))
        call check(error, network_parameter_count(2, 0), 3)
        if (allocated(error)) return
        call check(error, network_parameter_count(2, 2), 9)
        if (allocated(error)) return
        call network_predict(2, 0, [0.5_dp, -1.0_dp, 2.0_dp], &
            features, values, status)
        call check(error, status, profit_success)
        if (allocated(error)) return
        call check(error, maxval(abs(values - [3.5_dp, 0.25_dp, 0.25_dp])) < 1e-14_dp)
        if (allocated(error)) return
        theta = [0.3_dp, -0.2_dp, 0.4_dp, 0.5_dp, -0.1_dp, 0.2_dp, &
            0.7_dp, -0.6_dp, 0.3_dp]
        do i = 1, 3
            expected(i) = 0.3_dp + 0.7_dp*tanh( &
                0.3_dp*features(1, i) - 0.2_dp*features(2, i) - 0.1_dp) - &
                0.6_dp*tanh(0.4_dp*features(1, i) + &
                0.5_dp*features(2, i) + 0.2_dp)
        end do
        call network_predict(2, 2, theta, features, values, status)
        call check(error, status, profit_success)
        if (allocated(error)) return
        call check(error, maxval(abs(values - expected)) < 1e-14_dp)
        if (allocated(error)) return
        call network_predict(1, 1, [1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp], &
            reshape([1000.0_dp, -1000.0_dp, 0.0_dp], [1, 3]), &
            values, status)
        call check(error, status, profit_success)
        if (allocated(error)) return
        call check(error, all(values == [1.0_dp, -1.0_dp, 0.0_dp]))
    end subroutine test_values

    subroutine test_products(error)
        type(error_type), allocatable, intent(out) :: error
        real(dp), parameter :: steps(3) = [1e-4_dp, 1e-5_dp, 1e-6_dp]
        real(dp) :: features(2, 4), direction(2, 4), weight(4)
        real(dp) :: theta(9), theta_direction(9), theta_adjoint(9)
        real(dp) :: shifted_theta(9), shifted_features(2, 4)
        real(dp) :: feature_adjoint(2, 4), tangent(4), plus(4), minus(4)
        real(dp) :: difference, left, right, h
        integer :: i, j, status, nhidden, nparameter
        features = reshape([(sin(real(i, dp)), i=1, 8)], shape(features))
        direction = reshape([(cos(real(i, dp)), i=1, 8)], shape(direction))
        weight = [0.2_dp, -0.5_dp, 1.1_dp, -0.8_dp]
        theta = [(0.4_dp*sin(real(i, dp)), i=1, 9)]
        theta_direction = [(0.3_dp*cos(real(i, dp)), i=1, 9)]
        do nhidden = 0, 2, 2
            nparameter = network_parameter_count(2, nhidden)
            call network_jvp(2, nhidden, theta(:nparameter), features, &
                theta_direction(:nparameter), direction, tangent, status)
            call check(error, status, profit_success)
            if (allocated(error)) return
            do j = 1, size(steps)
                h = steps(j)
                shifted_theta = theta + h*theta_direction
                shifted_features = features + h*direction
                call network_predict(2, nhidden, shifted_theta(:nparameter), &
                    shifted_features, plus, status)
                call check(error, status, profit_success)
                if (allocated(error)) return
                shifted_theta = theta - h*theta_direction
                shifted_features = features - h*direction
                call network_predict(2, nhidden, shifted_theta(:nparameter), &
                    shifted_features, minus, status)
                call check(error, status, profit_success)
                if (allocated(error)) return
                difference = maxval(abs(tangent - (plus - minus)/(2*h)))
                call check(error, difference < 3e-9_dp)
                if (allocated(error)) return
            end do
            call network_vjp(2, nhidden, theta(:nparameter), features, weight, &
                theta_adjoint(:nparameter), feature_adjoint, status)
            call check(error, status, profit_success)
            if (allocated(error)) return
            left = dot_product(weight, tangent)
            right = dot_product(theta_direction(:nparameter), &
                theta_adjoint(:nparameter)) + &
                sum(direction*feature_adjoint)
            call check(error, left, right, thr=2e-14_dp)
            if (allocated(error)) return
            do i = 1, nparameter
                h = 1e-5_dp
                theta(i) = theta(i) + h
                call network_predict(2, nhidden, theta(:nparameter), &
                    features, plus, status)
                theta(i) = theta(i) - 2*h
                call network_predict(2, nhidden, theta(:nparameter), &
                    features, minus, status)
                theta(i) = theta(i) + h
                call check(error, theta_adjoint(i), &
                    dot_product(weight, plus - minus)/(2*h), thr=3e-10_dp)
                if (allocated(error)) return
            end do
        end do
    end subroutine test_products

    subroutine test_covariance(error)
        type(error_type), allocatable, intent(out) :: error
        integer, parameter :: n = 100000
        real(dp), parameter :: epsilon = 1e-3_dp
        real(dp) :: theta(9), centre(2, 1), directions(2, 1), zero_theta(9)
        real(dp) :: tangent(1), jacobian(2), predicted, empirical, average
        real(dp), allocatable :: features(:, :), values(:), uniform(:, :)
        real(dp) :: radius, z1, z2
        integer, allocatable :: seed(:)
        integer :: nseed, i, status
        theta = [0.3_dp, -0.2_dp, 0.4_dp, 0.5_dp, -0.1_dp, 0.2_dp, &
            0.7_dp, -0.6_dp, 0.3_dp]
        centre(:, 1) = [0.4_dp, -0.3_dp]
        zero_theta = 0.0_dp
        do i = 1, 2
            directions = 0.0_dp
            directions(i, 1) = 1.0_dp
            call network_jvp(2, 2, theta, centre, zero_theta, directions, &
                tangent, status)
            call check(error, status, profit_success)
            if (allocated(error)) return
            jacobian(i) = tangent(1)
        end do
        predicted = epsilon**2*((jacobian(1) + 0.4_dp*jacobian(2))**2 + &
            (0.8_dp*jacobian(2))**2)
        allocate (features(2, n), values(n), uniform(2, n))
        call random_seed(size=nseed)
        allocate (seed(nseed))
        seed = 72319
        call random_seed(put=seed)
        call random_number(uniform)
        do i = 1, n
            radius = sqrt(-2*log(max(uniform(1, i), tiny(1.0_dp))))
            z1 = radius*cos(2*acos(-1.0_dp)*uniform(2, i))
            z2 = radius*sin(2*acos(-1.0_dp)*uniform(2, i))
            features(1, i) = centre(1, 1) + epsilon*z1
            features(2, i) = centre(2, 1) + epsilon*(0.4_dp*z1 + 0.8_dp*z2)
        end do
        call network_predict(2, 2, theta, features, values, status)
        call check(error, status, profit_success)
        if (allocated(error)) return
        average = sum(values)/n
        empirical = sum((values - average)**2)/(n - 1)
        call check(error, abs(empirical/predicted - 1.0_dp) < 0.025_dp)
    end subroutine test_covariance

    subroutine test_positive(error)
        type(error_type), allocatable, intent(out) :: error
        real(dp) :: scores(9), values(9), derivatives(9), plus(9), minus(9), unused(9)
        real(dp) :: shifted_scores(9)
        integer :: status
        scores = [-1000.0_dp, -100.0_dp, -40.0_dp, -2.0_dp, 0.0_dp, &
            2.0_dp, 40.0_dp, 100.0_dp, 1000.0_dp]
        call network_positive_map(scores, values, derivatives, status)
        call check(error, status, profit_success)
        if (allocated(error)) return
        call check(error, all(ieee_is_finite(values)) .and. all(values >= 0.0_dp))
        if (allocated(error)) return
        call check(error, values(3)/exp(-40.0_dp), 1.0_dp, thr=1e-14_dp)
        if (allocated(error)) return
        call check(error, derivatives(3)/exp(-40.0_dp), 1.0_dp, thr=1e-14_dp)
        if (allocated(error)) return
        call check(error, values(5), log(2.0_dp), thr=1e-14_dp)
        if (allocated(error)) return
        call check(error, derivatives(5), 0.5_dp, thr=1e-14_dp)
        if (allocated(error)) return
        shifted_scores = scores + 1e-4_dp
        call network_positive_map(shifted_scores, plus, unused, status)
        shifted_scores = scores - 1e-4_dp
        call network_positive_map(shifted_scores, minus, unused, status)
        call check(error, maxval(abs(derivatives - (plus - minus)/2e-4_dp)) < 1e-9_dp)
        if (allocated(error)) return
        call network_positive_map([huge(1.0_dp), -huge(1.0_dp)], &
            values(:2), derivatives(:2), status)
        call check(error, status, profit_success)
        if (allocated(error)) return
        call check(error, all(values(:2) == [huge(1.0_dp), 0.0_dp]))
    end subroutine test_positive

    subroutine test_invalid(error)
        type(error_type), allocatable, intent(out) :: error
        real(dp) :: features(2, 1), theta(9), values(1), tangent(1), &
            theta_adjoint(9), feature_adjoint(2, 1), derivatives(1)
        integer :: status
        features = 0.0_dp
        theta = 0.0_dp
        call check(error, network_parameter_count(0, 2), 0)
        if (allocated(error)) return
        call check(error, network_parameter_count(2, -1), 0)
        if (allocated(error)) return
        call check(error, network_parameter_count(huge(1), huge(1)), 0)
        if (allocated(error)) return
        call network_predict(2, 2, theta(:8), features, values, status)
        call check(error, status, profit_invalid_shape)
        if (allocated(error)) return
        call check(error, all(values == 0.0_dp))
        if (allocated(error)) return
        theta(1) = ieee_value(0.0_dp, ieee_quiet_nan)
        call network_predict(2, 2, theta, features, values, status)
        call check(error, status, profit_invalid_input)
        if (allocated(error)) return
        theta = 0.0_dp
        call network_jvp(2, 2, theta, features, theta(:8), features, tangent, status)
        call check(error, status, profit_invalid_shape)
        if (allocated(error)) return
        call network_vjp(2, 2, theta, features, values, theta_adjoint(:8), &
            feature_adjoint, status)
        call check(error, status, profit_invalid_shape)
        if (allocated(error)) return
        theta(1) = huge(1.0_dp)
        features(1, 1) = huge(1.0_dp)
        call network_predict(2, 2, theta, features, values, status)
        call check(error, status, profit_range_error)
        if (allocated(error)) return
        call check(error, all(values == 0.0_dp))
        if (allocated(error)) return
        values(1) = ieee_value(0.0_dp, ieee_quiet_nan)
        call network_positive_map(values, tangent, derivatives, status)
        call check(error, status, profit_invalid_input)
    end subroutine test_invalid

    subroutine test_parallel(error)
        type(error_type), allocatable, intent(out) :: error
        integer, parameter :: n = 20000
        real(dp), allocatable :: features(:, :), direction(:, :), values(:)
        real(dp), allocatable :: reference(:), adjoint(:, :)
        real(dp), allocatable :: adjoint_reference(:, :), weight(:)
        real(dp) :: theta(9), theta_direction(9), gradient(9), gradient_reference(9)
        integer :: i, status
        allocate (features(2, n), direction(2, n), values(n), reference(n), &
            adjoint(2, n), adjoint_reference(2, n), weight(n))
        do i = 1, n
            features(1, i) = sin(0.01_dp*i)
            features(2, i) = cos(0.03_dp*i)
            direction(1, i) = cos(0.02_dp*i)
            direction(2, i) = sin(0.03_dp*i)
            weight(i) = cos(0.01_dp*i)
        end do
        theta = [(0.4_dp*sin(real(i, dp)), i=1, 9)]
        theta_direction = [(0.3_dp*cos(real(i, dp)), i=1, 9)]
        call omp_set_num_threads(1)
        call network_predict(2, 2, theta, features, reference, status)
        call omp_set_num_threads(4)
        call network_predict(2, 2, theta, features, values, status)
        call check(error, status, profit_success)
        if (allocated(error)) return
        call check(error, all(values == reference))
        if (allocated(error)) return
        call omp_set_num_threads(1)
        call network_jvp(2, 2, theta, features, theta_direction, direction, &
            reference, status)
        call omp_set_num_threads(4)
        call network_jvp(2, 2, theta, features, theta_direction, direction, &
            values, status)
        call check(error, status, profit_success)
        if (allocated(error)) return
        call check(error, all(values == reference))
        if (allocated(error)) return
        call omp_set_num_threads(1)
        call network_vjp(2, 2, theta, features, weight, gradient_reference, &
            adjoint_reference, status)
        call omp_set_num_threads(4)
        call network_vjp(2, 2, theta, features, weight, gradient, adjoint, status)
        call check(error, status, profit_success)
        if (allocated(error)) return
        call check(error, all(adjoint == adjoint_reference))
        if (allocated(error)) return
        call check(error, maxval(abs(gradient - gradient_reference))/ &
            max(1.0_dp, maxval(abs(gradient_reference))) < 1e-12_dp)
    end subroutine test_parallel
end module network_tests

program test_network
    use, intrinsic :: iso_fortran_env, only: error_unit
    use omp_lib, only: omp_set_num_threads
    use testdrive, only: run_testsuite
    use network_tests, only: collect
    implicit none
    integer :: status
    status = 0
    call omp_set_num_threads(1)
    call run_testsuite(collect, error_unit, status)
    if (status /= 0) error stop 1
end program test_network
