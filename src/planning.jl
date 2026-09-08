"""One executable or future solver route known to the planner."""
struct RouteCapability
    route::Symbol
    family::Symbol
    implemented::Bool
end

"""Mathematical eligibility decision for one candidate route."""
struct EligibilityDecision
    route::Symbol
    eligible::Bool
    reason::Symbol
end

"""Interpretation result for one policy layer."""
struct LayerDecision
    layer::Symbol
    choice::RouteChoice
    accepted::Bool
    reason::Symbol
end

"""Pure planning result. It contains no factorization, allocation of work buffers, or solve result."""
struct RoutePlan
    candidate_routes::Vector{Symbol}
    eligibility::Vector{EligibilityDecision}
    planned_routes::Vector{Symbol}
    execution_routes::Vector{Symbol}
    unavailable_routes::Vector{Symbol}
    layer_decisions::Vector{LayerDecision}
end

"""Opt-in, compact audit record of route qualification, attempts, and residual acceptance."""
struct RouteCertificate
    contract::MathematicalContract
    candidate_routes::Vector{Symbol}
    qualified_routes::Vector{Symbol}
    unavailable_routes::Vector{Symbol}
    layer_decisions::Vector{LayerDecision}
    attempted_routes::Vector{Symbol}
    selected_route::Union{Nothing, Symbol}
    fallback_reason::Union{Nothing, Symbol}
    residual_norm::Union{Nothing, Float64}
    residual_ratio::Union{Nothing, Float64}
    residual_accepted::Bool
    notes::Vector{String}
end

const _DIRECT_CAPABILITIES = (
    RouteCapability(:cholesky, :direct, true),
    RouteCapability(:lu, :direct, true),
    RouteCapability(:qr, :direct, true),
    RouteCapability(:svd, :direct, true),
)

const _ITERATIVE_CAPABILITIES = (
    RouteCapability(:cg, :iterative, true),
    RouteCapability(:minres, :iterative, true),
    RouteCapability(:gmres, :iterative, true),
    RouteCapability(:fgmres, :iterative, true),
    RouteCapability(:bicgstab, :iterative, true),
)

_routes(capabilities) = Symbol[capability.route for capability in capabilities]
_is_square(A) = size(A, 1) == size(A, 2)

function _choice_decision(layer::Symbol, choice::RouteChoice, supported::Set{Symbol})
    choice isa Auto && return LayerDecision(layer, choice, true, :automatic)
    choice isa Forbid && return LayerDecision(layer, choice, true, :forbidden)
    choice.route in supported && return LayerDecision(layer, choice, true, :supported)
    reason = choice isa Lock ? :locked_choice_not_implemented : :preferred_choice_not_implemented
    return LayerDecision(layer, choice, false, reason)
end

function _allows(choice::RouteChoice, value::Symbol)
    choice isa Auto && return true
    choice isa Prefer && return true
    choice isa Lock && return choice.route == value
    choice isa Forbid && return choice.route != value
    return false
end

function _prefer_first(routes::Vector{Symbol}, choice::RouteChoice)
    choice isa Prefer || return routes
    choice.route in routes || return routes
    return vcat(Symbol[choice.route], [route for route in routes if route != choice.route])
end

function _direct_order(problem::AdaptiveLinearProblem)
    if _certified_or_proved(problem.contract.rank_deficient)
        return Symbol[:svd, :qr, :lu, :cholesky]
    end
    return Symbol[:cholesky, :lu, :qr, :svd]
end

function _default_iterative_route(problem::AdaptiveLinearProblem, forbidden::Union{Nothing, Symbol}=nothing)
    candidates = if _certified_or_proved(problem.contract.hermitian) &&
                    _certified_or_proved(problem.contract.positive_definite)
        Symbol[:cg, :minres, :gmres, :fgmres]
    elseif _certified_or_proved(problem.contract.hermitian)
        Symbol[:minres, :gmres, :fgmres]
    elseif _fixed_linear_preconditioner(problem)
        Symbol[:gmres, :fgmres]
    else
        Symbol[:fgmres]
    end
    for method in candidates
        method == forbidden && continue
        qualify_iterative(problem, method).eligible && return method
    end
    return nothing
end

function _iterative_candidates(problem::AdaptiveLinearProblem, policy::RoutePolicy)
    choice = policy.iterative
    if choice isa Lock || choice isa Prefer
        routes = Symbol[choice.route]
        if choice isa Prefer && choice.route == :gmres && !_fixed_linear_preconditioner(problem)
            push!(routes, :fgmres)
        end
        return routes
    elseif choice isa Forbid
        method = _default_iterative_route(problem, choice.route)
        return method === nothing ? Symbol[] : Symbol[method]
    end
    method = _default_iterative_route(problem)
    return method === nothing ? Symbol[] : Symbol[method]
end

function _eligibility(route::Symbol, problem::AdaptiveLinearProblem)
    contract = problem.contract
    problem.A isa AbstractMatrix ||
        return EligibilityDecision(route, false, :matrix_free_direct_route_unavailable)
    if route == :cholesky
        _is_square(problem.A) || return EligibilityDecision(route, false, :nonsquare)
        _certified_or_proved(contract.hermitian) || return EligibilityDecision(route, false, :hermitian_evidence_insufficient)
        _certified_or_proved(contract.positive_definite) || return EligibilityDecision(route, false, :positive_definite_evidence_insufficient)
        return EligibilityDecision(route, true, :qualified)
    elseif route == :lu
        _is_square(problem.A) || return EligibilityDecision(route, false, :nonsquare)
        _certified_or_proved(contract.rank_deficient) && return EligibilityDecision(route, false, :known_rank_deficient)
        return EligibilityDecision(route, true, :qualified)
    elseif route in (:qr, :svd)
        return EligibilityDecision(route, true, :qualified)
    end
    return EligibilityDecision(route, false, :unknown_route)
end

"""
    plan(problem, policy=RoutePolicy())

Build a route plan without executing numerical kernels. In version 0.0.2 only the
direct family is executable; requests that lock an unavailable iterative or
preconditioner layer produce an empty executable route list and an explicit decision.
"""
function plan(problem::AdaptiveLinearProblem, policy::RoutePolicy=RoutePolicy())
    direct_routes = _routes(_DIRECT_CAPABILITIES)
    iterative_routes = _routes(_ITERATIVE_CAPABILITIES)
    layer_decisions = LayerDecision[
        _choice_decision(:family, policy.family, Set((:direct, :iterative))),
        _choice_decision(:direct, policy.direct, Set(direct_routes)),
        _choice_decision(:iterative, policy.iterative, Set(iterative_routes)),
        _choice_decision(:preconditioner, policy.preconditioner, Set{Symbol}()),
        _choice_decision(:fallback, policy.fallback, Set{Symbol}()),
    ]

    direct_enabled = _allows(policy.family, :direct) && !(policy.iterative isa Lock)
    direct_candidates = direct_enabled ? [route for route in _direct_order(problem)
        if _allows(policy.direct, route)] : Symbol[]
    direct_candidates = _prefer_first(direct_candidates, policy.direct)

    direct_eligibility = [_eligibility(route, problem) for route in direct_candidates]
    direct_route_available = any(decision -> decision.eligible, direct_eligibility)

    iterative_requested = policy.family isa Lock || policy.family isa Prefer ||
                          policy.iterative isa Lock || policy.iterative isa Prefer ||
                          (!direct_route_available && !(policy.direct isa Lock) &&
                           !(policy.family isa Lock && policy.family.route == :direct))
    iterative_enabled = iterative_requested && _allows(policy.family, :iterative)
    requested_iterative = iterative_enabled ? _iterative_candidates(problem, policy) : Symbol[]
    iterative_candidates = [route for route in requested_iterative if _allows(policy.iterative, route)]

    candidates = if policy.family isa Lock && policy.family.route == :iterative
        iterative_candidates
    elseif policy.iterative isa Lock || policy.iterative isa Prefer || policy.family isa Prefer
        vcat(iterative_candidates, direct_candidates)
    else
        vcat(direct_candidates, iterative_candidates)
    end
    eligibility = EligibilityDecision[]
    for route in candidates
        if route in direct_routes
            push!(eligibility, _eligibility(route, problem))
        else
            qualification = qualify_iterative(problem, route)
            push!(eligibility, EligibilityDecision(
                qualification.method, qualification.eligible, qualification.reason))
        end
    end
    planned_routes = Symbol[decision.route for decision in eligibility if decision.eligible]
    iterative_backend_available = _iterative_backend_available(problem)
    execution_routes = Symbol[route for route in planned_routes if
        route in direct_routes || (route in iterative_routes && iterative_backend_available)]
    unavailable_routes = Symbol[route for route in planned_routes if !(route in execution_routes)]

    preconditioner_locked = policy.preconditioner isa Lock
    fallback_locked = policy.fallback isa Lock
    if preconditioner_locked || fallback_locked
        execution_routes = Symbol[]
    end

    return RoutePlan(candidates, eligibility, planned_routes, execution_routes, unavailable_routes, layer_decisions)
end
