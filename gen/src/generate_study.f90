program generate_study
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, num, operator(+), operator(-), &
        operator(*), operator(/), operator(**), operator(==), &
        exp, log, sqrt, erfc, pi_expr
    use fortsym_string, only: str, chars
    use fortsym_kernel, only: kernel_spec_t, emit_kernel
    use fortsym_diff, only: diff
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    implicit none
    type(arena_t), target :: arena
    type(expr_t) :: residual, sigma, z, pi, nll, improve, probability, density
    type(expr_t) :: value, roots(3), scale, squared, count, critical, average, se
    type(symengine_engine_t) :: engine
    type(engine_result_t) :: simplified
    type(kernel_spec_t) :: spec
    character(4096) :: directory
    integer :: i

    if (command_argument_count() /= 1) error stop "supply output directory"
    call get_command_argument(1, directory)
    call arena%init()
    engine = make_symengine_engine(arena)
    pi = pi_expr(arena)
    residual = sym(arena, "residual")
    sigma = sym(arena, "sigma")
    z = residual/sigma
    nll = log(sigma) + log(2*pi)/2 + z**2/2
    call configure("score_leaf", [str("residual"), str("sigma")], &
        [str("standardized"), str("nll")])
    call emit([z, nll])

    improve = sym(arena, "improvement")
    z = improve/sigma
    probability = erfc(-z/sqrt(num(arena, 2)))/2
    density = exp(-z**2/2)/sqrt(2*pi)
    call configure("ei_leaf", [str("improvement"), str("sigma")], &
        [str("value"), str("d_improvement"), str("d_sigma")])
    value = improve*probability + sigma*density
    roots(1:2) = [diff(value, improve) - probability, diff(value, sigma) - density]
    do i = 1, 2
        simplified = engine%simplify(roots(i))
        if (.not. simplified%ok) error stop "EI simplification failed"
        if (.not. (simplified%value == num(arena, 0))) &
            error stop "EI derivative identity failed"
    end do
    ! Preserve erfc(-z) and the standardized square in the emitted evaluation.
    call emit([value, probability, density])

    scale = sym(arena, "scale")
    squared = sym(arena, "squared")
    count = sym(arena, "count")
    critical = sym(arena, "critical")
        average = sym(arena, "average")
        call configure("square_leaf", [str("residual"), str("scale")], [str("value")])
        call emit([(residual/scale)**2])
        call configure("rmse_leaf", [str("scale"), str("squared"), str("count")], &
            [str("value")])
        call emit([scale*sqrt(squared/count)])
        se = scale*sqrt(squared/(count - 1)/count)
        call configure("replicate_leaf", [str("average"), str("scale"), &
            str("squared"), str("count"), str("critical")], &
            [str("se"), str("lower"), str("upper")])
        call emit([se, average - critical*se, average + critical*se])

    contains

        subroutine configure(name, args, outputs)
            use fortsym_string, only: str_t
            character(*), intent(in) :: name
            type(str_t), intent(in) :: args(:), outputs(:)
            type(kernel_spec_t) :: fresh
            spec = fresh
            spec%name = str(name)
            spec%module_name = str("profit_generated_"//name)
            spec%args = args
            spec%outputs = outputs
            spec%generator = str("gen/src/generate_study.f90")
            spec%generator_revision = str("063008070086aaf2d78f6e9d9a8fddff44c74a46")
            spec%regenerate_command = str("cd gen && fpm run -- ../src/generated")
            spec%elemental_procedure = .true.
        end subroutine configure

        subroutine emit(roots)
            type(expr_t), intent(in) :: roots(:)
            integer :: unit, ios
            open (newunit=unit, file=trim(directory)//"/"// &
                chars(spec%module_name)//".f90", status="replace", &
                access="stream", form="unformatted", iostat=ios)
            if (ios /= 0) error stop "cannot open generated output"
            write (unit, iostat=ios) chars(emit_kernel(roots, spec))
            if (ios /= 0) error stop "cannot write generated output"
            close (unit)
        end subroutine emit
    end program generate_study
