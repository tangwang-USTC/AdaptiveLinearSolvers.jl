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
        locked in _SUPPORTED_ROUTES || throw(ArgumentError("route $locked is not implemented in v0.1.1"))
        locked in forbidden && throw(ArgumentError("route $locked is both locked and forbidden"))
        _qualified(locked, problem) || throw(ArgumentError("locked route $locked is not mathematically qualified"))
        return Symbol[locked]
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
        route in _SUPPORTED_ROUTES || throw(ArgumentError("preferred route $route is not implemented in v0.1.1"))
        _qualified(route, problem) && !(route in default) && pushfirst!(default, route)
    end
    return [route for route in unique(default) if !(route in forbidden)]
end

function _solve_route(route::Symbol, A, b)
    route == :direct && return A \ b
    route == :lu && return lu(A) \ b
    route == :cholesky && return cholesky(Hermitian(A)) \ b
    route == :qr && return qr(A) \ b
    route == :svd && return svd(A) \ b
    throw(ArgumentError("route $route is not implemented in v0.1.1"))
end

function _residual_ratio(A, x, b)
    denominator = norm(b)
    numerator = norm(b - A * x)
    return denominator == 0 ? Float64(numerator) : Float64(numerator / denominator)
end

function _telemetry(problem::AdaptiveLinearProblem, route::Symbol, residual_ratio::Float64,
        telemetry::TelemetryPolicy, history::Union{Nothing, HistoryStore}, status::Symbol)
    telemetry.level == :off && return (nothing, nothing)
    basic = (route=telemetry.output.route ? route : nothing,
             residual_ratio=telemetry.output.residual_ratio ? residual_ratio : nothing,
             conditioning=telemetry.output.conditioning ? problem.conditioning : nothing)
    telemetry.level == :basic && return (basic, nothing)
    fingerprint = _fingerprint(problem)
    record = (fingerprint=fingerprint, route=route, residual_ratio=residual_ratio,
              status=status, label=problem.label)
    should_emit = isempty(telemetry.emit_on) || status in telemetry.emit_on
    should_emit && history !== nothing && _record!(history, record)
    return (merge(basic, (fingerprint=fingerprint, emitted=should_emit)), record)
end

"""
    solve(problem; policy=RoutePolicy(), telemetry=TelemetryPolicy(), history=nothing)

Solve an explicit dense or sparse linear system with a mathematically qualified direct route.
Version 0.1.1 deliberately does not infer symmetry or positive definiteness from samples.
"""
function solve(problem::AdaptiveLinearProblem;
        policy::RoutePolicy=RoutePolicy(), telemetry::TelemetryPolicy=TelemetryPolicy(),
        history::Union{Nothing, HistoryStore}=nothing)
    _validate_telemetry(telemetry)
    routes = _route_order(problem, policy)
    isempty(routes) && throw(ArgumentError("no permitted and qualified route remains"))
    last_error = nothing
    for route in routes
        try
            x = _solve_route(route, problem.A, problem.b)
            residual_ratio = _residual_ratio(problem.A, x, problem.b)
            telemetry_data, record = _telemetry(problem, route, residual_ratio, telemetry, history, :success)
            return AdaptiveLinearSolution(x, :success, route, residual_ratio, telemetry_data, record)
        catch error
            last_error = error
        end
    end
    throw(ErrorException("all qualified direct routes failed; last error: $(sprint(showerror, last_error))"))
end

solve(A, b; kwargs...) = solve(AdaptiveLinearProblem(A, b); kwargs...)
