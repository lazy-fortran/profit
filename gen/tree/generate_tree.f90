program generate_tree
    use fortsym_expr, only: expr_t, sym, num, operator(+), operator(-), &
        operator(*), operator(/), operator(**), operator(==)
    use fortsym_string, only: str, chars, words
    use fortsym_kernel, only: kernel_spec_t, emit_kernel
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    implicit none
    type(expr_t) :: weighted_sum, weight_sum, ridge
    type(expr_t) :: mean, stationary
    type(symengine_engine_t) :: engine
    type(engine_result_t) :: simplified
    type(kernel_spec_t) :: spec
    character(4096) :: directory
    integer :: unit, ios

    if (command_argument_count() /= 1) error stop "supply output directory"
    call get_command_argument(1, directory)
    engine = make_symengine_engine()
    weighted_sum = sym("weighted_sum")
    weight_sum = sym("weight_sum")
    ridge = sym("ridge")
    mean = weighted_sum/(weight_sum + ridge)
    stationary = 2*(weight_sum + ridge)*mean - 2*weighted_sum
    simplified = engine%simplify(stationary)
    if (.not. simplified%ok) error stop "tree leaf identity failed"
    if (.not. (simplified%value == num(0))) &
        error stop "tree leaf stationarity failed"
    spec%name = str("tree_leaf")
    spec%module_name = str("profit_generated_tree_leaf")
    spec%args = words("weighted_sum weight_sum ridge")
    spec%outputs = words("value score")
    spec%elemental_procedure = .true.
    spec%generator = str("gen/tree/generate_tree.f90")
    spec%generator_revision = str("b9ea95c5a0475b016b1e149875a9bec9d8562533")
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
