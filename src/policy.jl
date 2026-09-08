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

"""Per-layer controls: family, direct method, iterative method, preconditioner, and fallback."""
Base.@kwdef struct RoutePolicy
    family::RouteChoice = Auto()
    direct::RouteChoice = Auto()
    iterative::RouteChoice = Auto()
    preconditioner::RouteChoice = Auto()
    fallback::RouteChoice = Auto()
end
