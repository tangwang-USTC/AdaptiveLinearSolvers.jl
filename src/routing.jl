"""
    solve(A, b; kwargs...)

Public convenience entry point for an explicit linear system. Keyword arguments are
forwarded unchanged to the contract-aware `solve(problem; ...)` implementation below.
"""
solve(A, b; kwargs...) = solve(AdaptiveLinearProblem(A, b); kwargs...)

function _solve_route(route::Symbol, A, b)
    route == :direct && return A \ b
    route == :lu && return lu(A) \ b
    route == :cholesky && return cholesky(Hermitian(A)) \ b
    route == :qr && return qr(A) \ b
    route == :svd && return svd(A) \ b
    throw(ArgumentError("route $route is not implemented in v0.0.4"))
end

function _execute_route(route::Symbol, problem::AdaptiveLinearProblem,
        backend::BackendSelection, residual_policy::ResidualPolicy,
        iteration_control::IterationControl, resource_budget::ResourceBudget)
    if backend.name == :stdlib
        return _solve_route(route, problem.A, problem.b), nothing
    elseif backend.name == :krylov
        operator = CountingOperator(problem.A; limit=resource_budget.max_operator_applications)
        x, report = _solve_krylov(route, operator, problem.b, residual_policy,
            iteration_control, problem.preconditioner)
        report = IterationReport(report.method, report.iterations, report.converged,
            report.backend_status, report.elapsed_seconds, report.residual_history,
            operator.applications)
        return x, report
    elseif backend.name == :iterativesolvers
        isfinite(resource_budget.max_seconds) &&
            throw(ArgumentError("IterativeSolvers adapter does not support a hard time budget"))
        return _solve_iterativesolvers(route, problem.A, problem.b, residual_policy,
            iteration_control, problem.preconditioner)
    elseif backend.name == :linearsolve
        (isfinite(resource_budget.max_seconds) || resource_budget.max_operator_applications > 0) &&
            throw(ArgumentError("LinearSolve adapter does not support hard time or operator budgets"))
        return _solve_linearsolve(route, problem.A, problem.b, residual_policy,
            iteration_control, problem.preconditioner)
    end
    throw(ArgumentError("backend $(backend.name) cannot execute route $route"))
end

function _residual_metrics(A, x, b, policy::ResidualPolicy)
    residual_norm = Float64(norm(b - A * x))
    rhs_norm = Float64(norm(b))
    residual_ratio = rhs_norm == 0 ? nothing : residual_norm / rhs_norm
    accepted = residual_norm <= policy.absolute_tolerance ||
               (residual_ratio !== nothing && residual_ratio <= policy.relative_tolerance)
    return residual_norm, residual_ratio, accepted
end

function _validate_residual_policy(policy::ResidualPolicy)
    policy.absolute_tolerance >= 0 || throw(ArgumentError("absolute_tolerance must be nonnegative"))
    policy.relative_tolerance >= 0 || throw(ArgumentError("relative_tolerance must be nonnegative"))
    return policy
end

function _certificate(problem::AdaptiveLinearProblem, telemetry::TelemetryPolicy,
        route_plan::RoutePlan, attempted::Vector{Symbol},
        selected::Union{Nothing, Symbol}, fallback_reason::Union{Nothing, Symbol},
        residual_norm::Union{Nothing, Float64}, residual_ratio::Union{Nothing, Float64},
        accepted::Bool, diagnosis::Union{Nothing, NumericalDiagnosis}, notes::Vector{String})
    telemetry.level == :off && return nothing
    telemetry.output.certificate || return nothing
    return RouteCertificate(problem.contract, copy(route_plan.candidate_routes),
        copy(route_plan.planned_routes), copy(route_plan.unavailable_routes),
        copy(route_plan.layer_decisions), copy(attempted),
        selected, fallback_reason, residual_norm, residual_ratio, accepted, diagnosis, copy(notes))
end

function _telemetry(problem::AdaptiveLinearProblem, route::Symbol, residual_ratio::Union{Nothing, Float64},
        diagnosis::NumericalDiagnosis, iteration, telemetry::TelemetryPolicy,
        history::Union{Nothing, HistoryStore}, status::SolveStatus)
    telemetry.level == :off && return (nothing, nothing)
    basic = (route=telemetry.output.route ? route : nothing,
             residual_ratio=telemetry.output.residual_ratio ? residual_ratio : nothing,
             conditioning=telemetry.output.conditioning ? diagnosis.conditioning.information : nothing,
             diagnostics=telemetry.output.diagnostics ? diagnosis : nothing)
    telemetry.level == :basic && return (basic, nothing)
    fingerprint = _fingerprint(problem, diagnosis.conditioning, telemetry.fingerprint)
    event = _status_event(status)
    preconditioner = problem.preconditioner
    preconditioner_name = preconditioner === nothing ? :none : preconditioner.name
    record = (fingerprint=fingerprint, route=route, residual_ratio=residual_ratio,
              status=status, label=problem.label, preconditioner=preconditioner_name,
              diagnosis=diagnosis)
    should_emit = isempty(telemetry.emit_on) || event in telemetry.emit_on
    should_emit && history !== nothing && _record!(history, record)
    return (merge(basic, (fingerprint=fingerprint, emitted=should_emit,
        trace=_trace_payload(iteration, telemetry))), record)
end

"""
    solve(problem; policy=RoutePolicy(), residual_policy=ResidualPolicy(), iteration_control=IterationControl(), conditioning_policy=ConditioningPolicy(), telemetry=TelemetryPolicy(), history=nothing)

Solve an explicit dense or sparse linear system with a mathematically qualified direct route.
Version 0.0.4 deliberately does not infer symmetry or positive definiteness from samples.
"""
function solve(problem::AdaptiveLinearProblem;
        policy::RoutePolicy=RoutePolicy(), residual_policy::ResidualPolicy=ResidualPolicy(),
        iteration_control::IterationControl=IterationControl(),
        resource_budget::ResourceBudget=ResourceBudget(),
        backend_policy::BackendPolicy=BackendPolicy(),
        conditioning_policy::ConditioningPolicy=ConditioningPolicy(),
        telemetry::TelemetryPolicy=TelemetryPolicy(),
        history::Union{Nothing, HistoryStore}=nothing,
        history_policy::HistoryPolicy=HistoryPolicy())
    _validate_telemetry(telemetry)
    _validate_residual_policy(residual_policy)
    _validate_iteration_control(iteration_control)
    _validate_resource_budget(resource_budget)
    _validate_conditioning_policy(conditioning_policy)
    route_plan = plan(problem, policy; conditioning_policy=conditioning_policy,
        history=history, fingerprint_profile=telemetry.fingerprint,
        history_policy=history_policy)
    isempty(route_plan.execution_routes) && begin
        notes = ["no permitted, implemented, and mathematically qualified route remains"]
        append!(notes, ["$(decision.layer): $(decision.reason)" for decision in route_plan.layer_decisions if !decision.accepted])
        diagnosis = diagnose(problem; policy=conditioning_policy,
            status=QualificationRejected)
        certificate = _certificate(problem, telemetry, route_plan, Symbol[], nothing,
            :no_qualified_route, nothing, nothing, false, diagnosis, notes)
        return AdaptiveLinearSolution(nothing, QualificationRejected, nothing, nothing,
            certificate, nothing, nothing, nothing)
    end

    attempted = Symbol[]
    notes = String[]
    terminal_status = NumericalFailure
    last_iteration = nothing
    for route in route_plan.execution_routes
        push!(attempted, route)
        try
            resources = assess_resources(problem, route, resource_budget)
            if !resources.eligible
                terminal_status = BudgetTerminated
                push!(notes, "$route: $(resources.reason)")
                continue
            end
            backend = select_backend(problem, route, backend_policy, resource_budget)
            if !backend.available
                push!(notes, "$route: $(backend.reason)")
                continue
            end
            effective_control = _budgeted_iteration_control(iteration_control, resource_budget)
            x, iteration = _execute_route(route, problem, backend, residual_policy,
                effective_control, resource_budget)
            last_iteration = iteration
            residual_norm, residual_ratio, accepted = _residual_metrics(problem.A, x, problem.b, residual_policy)
            if iteration !== nothing && !iteration.converged
                terminal_status = _iteration_failure_status(iteration, effective_control)
                push!(notes, "$route: $(iteration.backend_status)")
                continue
            end
            accepted || throw(ErrorException("residual acceptance failed"))
            status = length(attempted) == 1 ? Success : FallbackSuccess
            fallback_reason = status == FallbackSuccess ? :prior_route_failed : nothing
            diagnosis = diagnose(problem; policy=conditioning_policy,
                iteration=iteration, status=status)
            certificate = _certificate(problem, telemetry, route_plan,
                attempted, route, fallback_reason, residual_norm, residual_ratio, true, diagnosis, notes)
            telemetry_data, record = _telemetry(problem, route, residual_ratio, diagnosis,
                iteration, telemetry, history, status)
            return AdaptiveLinearSolution(x, status, route, residual_ratio, certificate,
                iteration, telemetry_data, record)
        catch error
            error isa OperatorApplicationBudgetExceeded && (terminal_status = BudgetTerminated)
            push!(notes, "$route: $(sprint(showerror, error))")
        end
    end
    diagnosis = diagnose(problem; policy=conditioning_policy,
        iteration=last_iteration, status=terminal_status)
    certificate = _certificate(problem, telemetry, route_plan,
        attempted, nothing, :all_routes_failed, nothing, nothing, false, diagnosis, notes)
    return AdaptiveLinearSolution(nothing, terminal_status, nothing, nothing,
        certificate, nothing, nothing, nothing)
end
