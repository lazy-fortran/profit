program network_scaling
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use omp_lib, only: omp_get_wtime, omp_set_num_threads
    use profit_models_network, only: network_parameter_count, network_predict, &
        network_vjp
    use profit_status, only: profit_success
    implicit none
    integer, parameter :: p = 14, hidden = 8, n = 262144, repetitions = 12
    real(dp), allocatable :: theta(:), features(:, :), values(:), reference(:)
    real(dp), allocatable :: gradient(:), reference_gradient(:)
    real(dp), allocatable :: input_adjoint(:, :), weights(:)
    real(dp) :: start, elapsed(2), forward(2), error
    integer :: i, j, r, mode, status, nparameter
    nparameter = network_parameter_count(p, hidden)
    allocate (theta(nparameter), features(p, n), &
        values(n), reference(n), gradient(nparameter), &
        reference_gradient(nparameter), &
        input_adjoint(p, n), weights(n))
    do i = 1, size(theta)
        theta(i) = sin(real(i, dp))*0.1_dp
    end do
    do i = 1, n
        do j = 1, p
            features(j, i) = sin(0.001_dp*i*j)
        end do
        weights(i) = cos(0.003_dp*i)/n
    end do
    do mode = 1, 2
        call omp_set_num_threads(1 + 3*(mode - 1))
        call network_predict(p, hidden, theta, features, values, status)
        if (status /= profit_success) error stop "forward failed"
        start = omp_get_wtime()
        do r = 1, repetitions
            call network_predict(p, hidden, theta, features, values, status)
            if (status /= profit_success) error stop "forward failed"
        end do
        forward(mode) = (omp_get_wtime() - start)/repetitions
        start = omp_get_wtime()
        do r = 1, repetitions
            call network_vjp(p, hidden, theta, features, weights, gradient, &
                input_adjoint, status)
            if (status /= profit_success) error stop "reverse failed"
        end do
        elapsed(mode) = (omp_get_wtime() - start)/repetitions
        if (mode == 1) then
            reference = values
            reference_gradient = gradient
        end if
    end do
    error = maxval(abs(gradient - reference_gradient))
    if (any(values /= reference)) error stop "threaded forward differs"
    if (error > 1e-12_dp) error stop "threaded gradient differs"
    write (*, '(a)') "operation,one_thread_seconds,four_thread_seconds,speedup"
    write (*, '(a,3(es16.8,:,","))') "predict,", forward, forward(1)/forward(2)
    write (*, '(a,3(es16.8,:,","))') "vjp,", elapsed, elapsed(1)/elapsed(2)
    write (*, '(a,es16.8)') "max_parameter_gradient_difference,", error
end program network_scaling
