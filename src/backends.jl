const _BACKEND_CAPABILITIES = (
    BackendCapability(:stdlib, :LinearAlgebra, Symbol[:cholesky, :lu, :qr, :svd],
        Symbol[:serial_cpu], true),
    BackendCapability(:krylov, :Krylov, Symbol[:cg, :minres, :gmres, :fgmres, :bicgstab],
        Symbol[:serial_cpu], true),
    BackendCapability(:iterativesolvers, :IterativeSolvers,
        Symbol[:cg, :gmres, :bicgstab], Symbol[:serial_cpu], false),
    BackendCapability(:linearsolve, :LinearSolve, Symbol[:cholesky, :lu, :qr, :svd, :cg, :gmres],
        Symbol[:serial_cpu, :gpu], false),
    BackendCapability(:petsc, :PETSc, Symbol[:cg, :minres, :gmres, :fgmres, :bicgstab],
        Symbol[:serial_cpu, :distributed_cpu], false),
    BackendCapability(:cuda, :CUDA, Symbol[:cg, :gmres, :bicgstab], Symbol[:gpu], false),
    BackendCapability(:mpi, :MPI, Symbol[:cg, :minres, :gmres, :fgmres, :bicgstab],
        Symbol[:distributed_cpu], false),
)

"""Return declared backend capabilities, including unavailable adapters for auditability."""
backend_capabilities() = collect(_BACKEND_CAPABILITIES)

function _validate_resource_budget(budget::ResourceBudget)
    budget.max_seconds > 0 || throw(ArgumentError("resource max_seconds must be positive"))
    budget.max_iterations >= 0 || throw(ArgumentError("resource max_iterations must be nonnegative"))
    budget.max_memory_bytes >= 0 || throw(ArgumentError("resource max_memory_bytes must be nonnegative"))
    budget.max_operator_applications >= 0 ||
        throw(ArgumentError("resource max_operator_applications must be nonnegative"))
    budget.execution in (:serial_cpu, :gpu, :distributed_cpu) ||
        throw(ArgumentError("resource execution must be :serial_cpu, :gpu, or :distributed_cpu"))
    return budget
end

function _minimum_workspace_bytes(problem::AdaptiveLinearProblem, route::Symbol)
    n = max(size(problem.A)...)
    scalar_bytes = try
        max(sizeof(eltype(problem.A)), 8)
    catch
        8
    end
    if route in (:cg, :minres, :gmres, :fgmres, :bicgstab)
        return 6 * n * scalar_bytes
    elseif problem.A isa StridedMatrix
        multiplier = route == :svd ? 4 : 2
        return multiplier * n * n * scalar_bytes
    end
    return nothing
end

"""Check a conservative minimum workspace bound before selecting an implemented backend."""
function assess_resources(problem::AdaptiveLinearProblem, route::Symbol,
        budget::ResourceBudget=ResourceBudget())
    _validate_resource_budget(budget)
    budget.execution == :serial_cpu ||
        return ResourceAssessment(false, :requested_execution_mode_unavailable, nothing)
    estimated = _minimum_workspace_bytes(problem, route)
    budget.max_memory_bytes == 0 && return ResourceAssessment(true, :unbounded_memory, estimated)
    estimated === nothing && return ResourceAssessment(false, :memory_budget_not_verifiable, nothing)
    estimated <= budget.max_memory_bytes ||
        return ResourceAssessment(false, :memory_budget_exceeded, estimated)
    return ResourceAssessment(true, :within_memory_budget, estimated)
end

function _capability(name::Symbol)
    return findfirst(capability -> capability.name == name, _BACKEND_CAPABILITIES)
end

"""Select only a registered and implemented backend that supports the route and execution mode."""
function select_backend(problem::AdaptiveLinearProblem, route::Symbol,
        policy::BackendPolicy=BackendPolicy(), budget::ResourceBudget=ResourceBudget())
    _validate_resource_budget(budget)
    candidates = if policy.requested == :auto
        _BACKEND_CAPABILITIES
    else
        index = _capability(policy.requested)
        index === nothing && return BackendSelection(nothing, false, :unknown_backend)
        (_BACKEND_CAPABILITIES[index],)
    end
    for capability in candidates
        route in capability.routes || continue
        budget.execution in capability.execution_modes || continue
        capability.implemented && return BackendSelection(capability.name, true, :selected)
        policy.requested == :auto && continue
        return BackendSelection(capability.name, false, :backend_adapter_not_implemented)
    end
    return BackendSelection(nothing, false, :no_implemented_backend_for_route)
end

function _budgeted_iteration_control(control::IterationControl, budget::ResourceBudget)
    max_iterations = if budget.max_iterations == 0
        control.max_iterations
    elseif control.max_iterations == 0
        budget.max_iterations
    else
        min(control.max_iterations, budget.max_iterations)
    end
    budget.max_operator_applications > 0 &&
        (max_iterations = max_iterations == 0 ? budget.max_operator_applications :
            min(max_iterations, budget.max_operator_applications))
    return IterationControl(max_iterations=max_iterations,
        max_seconds=min(control.max_seconds, budget.max_seconds), restart=control.restart,
        record_history=control.record_history)
end
