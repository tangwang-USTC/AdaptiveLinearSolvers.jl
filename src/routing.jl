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
    throw(ArgumentError("route $route is not implemented in v0.0.1"))
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
        accepted::Bool, notes::Vector{String})
    telemetry.level == :off && return nothing
    telemetry.output.certificate || return nothing
    return RouteCertificate(problem.contract, copy(route_plan.candidate_routes),
        copy(route_plan.planned_routes), copy(route_plan.unavailable_routes),
        copy(route_plan.layer_decisions), copy(attempted),
        selected, fallback_reason, residual_norm, residual_ratio, accepted, copy(notes))
end

function _telemetry(problem::AdaptiveLinearProblem, route::Symbol, residual_ratio::Union{Nothing, Float64},
        telemetry::TelemetryPolicy, history::Union{Nothing, HistoryStore}, status::SolveStatus)
    telemetry.level == :off && return (nothing, nothing)
    basic = (route=telemetry.output.route ? route : nothing,
             residual_ratio=telemetry.output.residual_ratio ? residual_ratio : nothing,
             conditioning=telemetry.output.conditioning ? problem.conditioning : nothing)
    telemetry.level == :basic && return (basic, nothing)
    fingerprint = _fingerprint(problem)
    event = _status_event(status)
    record = (fingerprint=fingerprint, route=route, residual_ratio=residual_ratio,
              status=status, label=problem.label)
    should_emit = isempty(telemetry.emit_on) || event in telemetry.emit_on
    should_emit && history !== nothing && _record!(history, record)
    return (merge(basic, (fingerprint=fingerprint, emitted=should_emit)), record)
end

"""
    solve(problem; policy=RoutePolicy(), residual_policy=ResidualPolicy(), telemetry=TelemetryPolicy(), history=nothing)

Solve an explicit dense or sparse linear system with a mathematically qualified direct route.
Version 0.0.1 deliberately does not infer symmetry or positive definiteness from samples.
"""
function solve(problem::AdaptiveLinearProblem;
        policy::RoutePolicy=RoutePolicy(), residual_policy::ResidualPolicy=ResidualPolicy(),
        telemetry::TelemetryPolicy=TelemetryPolicy(),
        history::Union{Nothing, HistoryStore}=nothing)
    _validate_telemetry(telemetry)
    _validate_residual_policy(residual_policy)
    route_plan = plan(problem, policy)
    isempty(route_plan.execution_routes) && begin
        notes = ["no permitted, implemented, and mathematically qualified route remains"]
        append!(notes, ["$(decision.layer): $(decision.reason)" for decision in route_plan.layer_decisions if !decision.accepted])
        certificate = _certificate(problem, telemetry, route_plan, Symbol[], nothing,
            :no_qualified_route, nothing, nothing, false, notes)
        return AdaptiveLinearSolution(nothing, QualificationRejected, nothing, nothing,
            certificate, nothing, nothing)
    end

    attempted = Symbol[]
    notes = String[]
    for route in route_plan.execution_routes
        push!(attempted, route)
        try
            x = _solve_route(route, problem.A, problem.b)
            residual_norm, residual_ratio, accepted = _residual_metrics(problem.A, x, problem.b, residual_policy)
            accepted || throw(ErrorException("residual acceptance failed"))
            status = length(attempted) == 1 ? Success : FallbackSuccess
            fallback_reason = status == FallbackSuccess ? :prior_route_failed : nothing
            certificate = _certificate(problem, telemetry, route_plan,
                attempted, route, fallback_reason, residual_norm, residual_ratio, true, notes)
            telemetry_data, record = _telemetry(problem, route, residual_ratio, telemetry, history, status)
            return AdaptiveLinearSolution(x, status, route, residual_ratio, certificate, telemetry_data, record)
        catch error
            push!(notes, "$route: $(sprint(showerror, error))")
        end
    end
    certificate = _certificate(problem, telemetry, route_plan,
        attempted, nothing, :all_routes_failed, nothing, nothing, false, notes)
    return AdaptiveLinearSolution(nothing, NumericalFailure, nothing, nothing,
        certificate, nothing, nothing)
end
