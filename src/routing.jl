_is_square(A) = size(A, 1) == size(A, 2)

function _qualified(route::Symbol, problem::AdaptiveLinearProblem)
    contract = problem.contract
    if route == :cholesky
        return _is_square(problem.A) && _certified_or_proved(contract.hermitian) &&
               _certified_or_proved(contract.positive_definite)
    elseif route in (:lu, :qr, :svd, :direct)
        return true
    end
    return false
end

function _route_order(problem::AdaptiveLinearProblem, policy::RoutePolicy)
    forbidden = _forbidden_routes(policy)
    locked = _locked_route(policy)
    if locked !== nothing
        locked in _SUPPORTED_ROUTES || throw(ArgumentError("route $locked is not implemented in v0.1.2"))
        locked in forbidden && throw(ArgumentError("route $locked is both locked and forbidden"))
        _qualified(locked, problem) || throw(ArgumentError("locked route $locked is not mathematically qualified"))
        return (candidate_routes=Symbol[locked], qualified_routes=Symbol[locked])
    end

    default = if _qualified(:cholesky, problem)
        Symbol[:cholesky, :lu, :qr, :svd]
    elseif _certified_or_proved(problem.contract.rank_deficient)
        Symbol[:svd, :qr, :lu]
    elseif _is_square(problem.A)
        Symbol[:lu, :qr, :svd]
    else
        Symbol[:qr, :svd]
    end
    for route in _preferred_routes(policy)
        route in _SUPPORTED_ROUTES || throw(ArgumentError("preferred route $route is not implemented in v0.1.2"))
        _qualified(route, problem) && !(route in default) && pushfirst!(default, route)
    end
    candidates = unique(default)
    return (candidate_routes=candidates,
            qualified_routes=[route for route in candidates if !(route in forbidden)])
end

function _solve_route(route::Symbol, A, b)
    route == :direct && return A \ b
    route == :lu && return lu(A) \ b
    route == :cholesky && return cholesky(Hermitian(A)) \ b
    route == :qr && return qr(A) \ b
    route == :svd && return svd(A) \ b
    throw(ArgumentError("route $route is not implemented in v0.1.2"))
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
        candidates::Vector{Symbol}, qualified::Vector{Symbol}, attempted::Vector{Symbol},
        selected::Union{Nothing, Symbol}, fallback_reason::Union{Nothing, Symbol},
        residual_norm::Union{Nothing, Float64}, residual_ratio::Union{Nothing, Float64},
        accepted::Bool, notes::Vector{String})
    telemetry.level == :off && return nothing
    telemetry.output.certificate || return nothing
    return RouteCertificate(problem.contract, copy(candidates), copy(qualified), copy(attempted),
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
Version 0.1.2 deliberately does not infer symmetry or positive definiteness from samples.
"""
function solve(problem::AdaptiveLinearProblem;
        policy::RoutePolicy=RoutePolicy(), residual_policy::ResidualPolicy=ResidualPolicy(),
        telemetry::TelemetryPolicy=TelemetryPolicy(),
        history::Union{Nothing, HistoryStore}=nothing)
    _validate_telemetry(telemetry)
    _validate_residual_policy(residual_policy)
    try
        plan = _route_order(problem, policy)
    catch error
        notes = [sprint(showerror, error)]
        certificate = _certificate(problem, telemetry, Symbol[], Symbol[], Symbol[], nothing,
            :qualification_rejected, nothing, nothing, false, notes)
        return AdaptiveLinearSolution(nothing, QualificationRejected, nothing, nothing,
            certificate, nothing, nothing)
    end
    isempty(plan.qualified_routes) && begin
        notes = ["no permitted and qualified route remains"]
        certificate = _certificate(problem, telemetry, plan.candidate_routes, plan.qualified_routes,
            Symbol[], nothing, :no_qualified_route, nothing, nothing, false, notes)
        return AdaptiveLinearSolution(nothing, QualificationRejected, nothing, nothing,
            certificate, nothing, nothing)
    end

    attempted = Symbol[]
    notes = String[]
    for route in plan.qualified_routes
        push!(attempted, route)
        try
            x = _solve_route(route, problem.A, problem.b)
            residual_norm, residual_ratio, accepted = _residual_metrics(problem.A, x, problem.b, residual_policy)
            accepted || throw(ErrorException("residual acceptance failed"))
            status = length(attempted) == 1 ? Success : FallbackSuccess
            fallback_reason = status == FallbackSuccess ? :prior_route_failed : nothing
            certificate = _certificate(problem, telemetry, plan.candidate_routes, plan.qualified_routes,
                attempted, route, fallback_reason, residual_norm, residual_ratio, true, notes)
            telemetry_data, record = _telemetry(problem, route, residual_ratio, telemetry, history, status)
            return AdaptiveLinearSolution(x, status, route, residual_ratio, certificate, telemetry_data, record)
        catch error
            push!(notes, "$route: $(sprint(showerror, error))")
        end
    end
    certificate = _certificate(problem, telemetry, plan.candidate_routes, plan.qualified_routes,
        attempted, nothing, :all_routes_failed, nothing, nothing, false, notes)
    return AdaptiveLinearSolution(nothing, NumericalFailure, nothing, nothing,
        certificate, nothing, nothing)
end

solve(A, b; kwargs...) = solve(AdaptiveLinearProblem(A, b); kwargs...)
