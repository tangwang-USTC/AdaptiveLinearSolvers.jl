abstract type RouteChoice end

struct Auto <: RouteChoice end
struct Prefer <: RouteChoice
    route::Symbol
end
struct Lock <: RouteChoice
    route::Symbol
end
struct Forbid <: RouteChoice
    route::Symbol
end

"""Per-layer controls. A single layer can be locked while all other layers remain automatic."""
Base.@kwdef struct RoutePolicy
    family::RouteChoice = Auto()
    direct::RouteChoice = Auto()
    iterative::RouteChoice = Auto()
    preconditioner::RouteChoice = Auto()
    fallback::RouteChoice = Auto()
end

const _SUPPORTED_ROUTES = Set((:direct, :lu, :cholesky, :qr, :svd))

function _forbidden_routes(policy::RoutePolicy)
    routes = Set{Symbol}()
    for choice in (policy.family, policy.direct, policy.iterative, policy.preconditioner, policy.fallback)
        choice isa Forbid && push!(routes, choice.route)
    end
    return routes
end

function _locked_route(policy::RoutePolicy)
    locks = Symbol[]
    for choice in (policy.family, policy.direct, policy.iterative, policy.preconditioner, policy.fallback)
        choice isa Lock && push!(locks, choice.route)
    end
    length(unique(locks)) <= 1 || throw(ArgumentError("conflicting locked routes: $(unique(locks))"))
    return isempty(locks) ? nothing : only(locks)
end

function _preferred_routes(policy::RoutePolicy)
    routes = Symbol[]
    for choice in (policy.family, policy.direct, policy.iterative, policy.preconditioner, policy.fallback)
        choice isa Prefer && push!(routes, choice.route)
    end
    return unique(routes)
end
