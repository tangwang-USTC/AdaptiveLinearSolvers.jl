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
    execution_routes::Vector{Symbol}
    layer_decisions::Vector{LayerDecision}
end

"""Opt-in, compact audit record of route qualification, attempts, and residual acceptance."""
struct RouteCertificate
    contract::MathematicalContract
    candidate_routes::Vector{Symbol}
    qualified_routes::Vector{Symbol}
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

function _eligibility(route::Symbol, problem::AdaptiveLinearProblem)
    contract = problem.contract
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

Build a route plan without executing numerical kernels. In version 0.0.1 only the
direct family is executable; requests that lock an unavailable iterative or
preconditioner layer produce an empty executable route list and an explicit decision.
"""
function plan(problem::AdaptiveLinearProblem, policy::RoutePolicy=RoutePolicy())
    direct_routes = _routes(_DIRECT_CAPABILITIES)
    layer_decisions = LayerDecision[
        _choice_decision(:family, policy.family, Set((:direct,))),
        _choice_decision(:direct, policy.direct, Set(direct_routes)),
        _choice_decision(:iterative, policy.iterative, Set{Symbol}()),
        _choice_decision(:preconditioner, policy.preconditioner, Set{Symbol}()),
        _choice_decision(:fallback, policy.fallback, Set{Symbol}()),
    ]

    candidates = [route for route in _direct_order(problem)
                  if _allows(policy.family, :direct) && _allows(policy.direct, route)]
    candidates = _prefer_first(candidates, policy.direct)
    eligibility = [_eligibility(route, problem) for route in candidates]
    execution_routes = Symbol[decision.route for decision in eligibility if decision.eligible]

    iterative_locked = policy.iterative isa Lock
    preconditioner_locked = policy.preconditioner isa Lock
    fallback_locked = policy.fallback isa Lock
    if iterative_locked || preconditioner_locked || fallback_locked
        execution_routes = Symbol[]
    end

    return RoutePlan(candidates, eligibility, execution_routes, layer_decisions)
end
