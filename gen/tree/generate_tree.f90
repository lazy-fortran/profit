program generate_tree
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, num, operator(+), operator(-), &
        operator(*), operator(/), operator(**), operator(==)
    use fortsym_string, only: str, chars
    use fortsym_kernel, only: kernel_spec_t, emit_kernel
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    implicit none
    type(arena_t), target :: arena
    type(expr_t) :: weighted_sum, weight_sum, ridge
    type(expr_t) :: mean, stationary
    type(symengine_engine_t) :: engine
    type(engine_result_t) :: simplified
    type(kernel_spec_t) :: spec
    character(4096) :: directory
    integer :: unit, ios

    if (command_argument_count() /= 1) error stop "supply output directory"
    call get_command_argument(1, directory)
    call arena%init()
    engine = make_symengine_engine(arena)
    weighted_sum = sym(arena, "weighted_sum")
    weight_sum = sym(arena, "weight_sum")
    ridge = sym(arena, "ridge")
    mean = weighted_sum/(weight_sum + ridge)
    stationary = 2*(weight_sum + ridge)*mean - 2*weighted_sum
    simplified = engine%simplify(stationary)
    if (.not. simplified%ok) error stop "tree leaf identity failed"
    if (.not. (simplified%value == num(arena, 0))) &
        error stop "tree leaf stationarity failed"
    spec%name = str("tree_leaf")
    spec%module_name = str("profit_generated_tree_leaf")
    spec%args = [str("weighted_sum"), str("weight_sum"), str("ridge")]
    spec%outputs = [str("value"), str("score")]
    spec%elemental_procedure = .true.
    spec%generator = str("gen/tree/generate_tree.f90")
    spec%generator_revision = str("063008070086aaf2d78f6e9d9a8fddff44c74a46")
    spec%regenerate_command = str("FORTSYM_BUILD=<built-tree> sh gen/tree/regenerate.sh")
    open(newunit=unit, file=trim(directory)//"/profit_generated_tree_leaf.f90", &
        status="replace", access="stream", form="unformatted", iostat=ios)
    if (ios /= 0) error stop "cannot open generated tree leaf"
    write(unit, iostat=ios) chars(emit_kernel([ &
        mean, &
        weighted_sum*weighted_sum/(weight_sum + ridge)], spec))
    if (ios /= 0) error stop "cannot write generated tree leaf"
    close(unit)
end program generate_tree
