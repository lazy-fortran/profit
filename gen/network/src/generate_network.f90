program generate_network
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, num, operator(+), operator(-), &
        operator(*), operator(/), operator(**), operator(==), &
        exp, log, tanh, abs
    use fortsym_string, only: str, chars, str_t
    use fortsym_kernel, only: kernel_spec_t, emit_kernel
    use fortsym_diff, only: diff
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    implicit none
    type(arena_t), target :: arena
    type(expr_t) :: z, t, e, rational_activation, slope, series
    type(symengine_engine_t) :: engine
    type(engine_result_t) :: identity
    type(kernel_spec_t) :: spec
    character(4096) :: directory
    integer :: j

    if (command_argument_count() /= 1) error stop "supply output directory"
    call get_command_argument(1, directory)
    call arena%init()
    engine = make_symengine_engine(arena)
    z = sym(arena, "z")
    t = sym(arena, "t")
    e = exp(-z)
    rational_activation = (1 - e**2)/(1 + e**2)
    slope = 4*e**2/(1 + e**2)**2
    identity = engine%simplify(diff(rational_activation, z) - slope)
    if (.not. identity%ok) error stop "activation derivative check failed"
    if (.not. (identity%value == num(arena, 0))) &
        error stop "activation derivative identity failed"
    e = exp(-abs(z))
    call configure("network_activation", [str("z")], &
        [str("value"), str("derivative")])
    call emit([tanh(z), 4*e**2/(1 + e**2)**2])

    call configure("network_positive", [str("z")], &
        [str("tail"), str("correction"), str("slope")])
    call emit([exp(-z), log(1 + exp(-z)), 1/(1 + exp(-z))])
    series = num(arena, 0)
    do j = 1, 7
        series = series + ((-1)**(j + 1))*t**j/num(arena, j)
    end do
    call configure("network_small_log", [str("t")], [str("value")])
    call emit([series])
contains
    subroutine configure(name, args, outputs)
        character(*), intent(in) :: name
        type(str_t), intent(in) :: args(:), outputs(:)
        type(kernel_spec_t) :: fresh
        spec = fresh
        spec%name = str(name)
        spec%module_name = str("profit_generated_"//name)
        spec%args = args
        spec%outputs = outputs
        spec%generator = str("gen/network/src/generate_network.f90")
        spec%generator_revision = str("063008070086aaf2d78f6e9d9a8fddff44c74a46")
        spec%regenerate_command = str("sh gen/network_regenerate.sh")
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
end program generate_network
