"""Request static fields from the current call; this does not create a fingerprint."""
Base.@kwdef struct OutputRequest
    route::Bool = false
    residual_ratio::Bool = false
    conditioning::Bool = false
    diagnostics::Bool = false
    certificate::Bool = false
end

"""Select labels needed for history-driven routing advice."""
Base.@kwdef struct FingerprintProfile
    representation::Bool = true
    size_band::Bool = true
    structure::Bool = true
    conditioning::Bool = true
    execution::Bool = true
end

"""Telemetry is off by default; expensive trace capture is bounded by `TelemetryBudget`."""
Base.@kwdef struct TelemetryPolicy
    level::Symbol = :off
    output::OutputRequest = OutputRequest()
    fingerprint::FingerprintProfile = FingerprintProfile()
    emit_on::Tuple{Vararg{Symbol}} = ()
    sample_every::Int = 1
    budget::TelemetryBudget = TelemetryBudget()
end

function _validate_telemetry(policy::TelemetryPolicy)
    policy.level in (:off, :basic, :fingerprint, :trace, :diagnostic) ||
        throw(ArgumentError("telemetry level $(policy.level) must be :off, :basic, :fingerprint, :trace, or :diagnostic"))
    policy.sample_every > 0 || throw(ArgumentError("telemetry sample_every must be positive"))
    policy.budget.max_trace_samples >= 0 ||
        throw(ArgumentError("telemetry max_trace_samples must be nonnegative"))
    policy.budget.max_extra_operator_applications >= 0 ||
        throw(ArgumentError("telemetry max_extra_operator_applications must be nonnegative"))
    policy.budget.max_seconds >= 0 ||
        throw(ArgumentError("telemetry max_seconds must be nonnegative"))
    return policy
end

Base.@kwdef struct MatrixFingerprint
    representation::Union{Nothing, Symbol} = nothing
    size_band::Union{Nothing, Symbol} = nothing
    structure::Union{Nothing, Symbol} = nothing
    conditioning::Union{Nothing, Symbol} = nothing
    execution::Union{Nothing, Symbol} = nothing
end

function _size_band(A)
    n = max(size(A)...)
    n <= 64 && return :small
    n <= 2_048 && return :medium
    return :large
end

function _structure_tag(contract::MathematicalContract)
    if _certified_or_proved(contract.hermitian) && _certified_or_proved(contract.positive_definite)
        return :hermitian_positive_definite
    elseif _certified_or_proved(contract.hermitian)
        return :hermitian
    elseif _certified_or_proved(contract.rank_deficient)
        return :rank_deficient
    end
    return :general
end

function _conditioning_tag(assessment::ConditioningAssessment)
    return assessment.state
end

function _fingerprint(problem::AdaptiveLinearProblem,
        assessment::ConditioningAssessment=assess_conditioning(problem),
        profile::FingerprintProfile=FingerprintProfile())
    return MatrixFingerprint(
        representation=profile.representation ? (problem.A isa AbstractLinearOperator ? :matrix_free :
            (issparse(problem.A) ? :sparse_explicit : :dense_explicit)) : nothing,
        size_band=profile.size_band ? _size_band(problem.A) : nothing,
        structure=profile.structure ? _structure_tag(problem.contract) : nothing,
        conditioning=profile.conditioning ? _conditioning_tag(assessment) : nothing,
        execution=profile.execution ? :serial_cpu : nothing,
    )
end

"""Build a low-cost, profile-controlled fingerprint without numerical estimation."""
function fingerprint(problem::AdaptiveLinearProblem;
        profile::FingerprintProfile=FingerprintProfile(),
        conditioning_policy::ConditioningPolicy=ConditioningPolicy())
    return _fingerprint(problem, assess_conditioning(problem, conditioning_policy), profile)
end

"""Bounded, opt-in in-memory history for fingerprint telemetry."""
mutable struct HistoryStore
    capacity::Int
    records::Vector{NamedTuple}
end

HistoryStore(capacity::Integer=256) = HistoryStore(Int(capacity), NamedTuple[])

function _record!(store::HistoryStore, record::NamedTuple)
    store.capacity > 0 || return store
    push!(store.records, record)
    while length(store.records) > store.capacity
        popfirst!(store.records)
    end
    return store
end

function _fingerprints_match(query::MatrixFingerprint, candidate::MatrixFingerprint)
    for field in fieldnames(MatrixFingerprint)
        query_value = getfield(query, field)
        query_value === nothing && continue
        getfield(candidate, field) == query_value || return false
    end
    return true
end

"""Return bounded in-memory records matching all fields enabled in `fingerprint`."""
function similar_records(store::HistoryStore, query::MatrixFingerprint;
        label::Union{Nothing, Symbol}=nothing)
    return [record for record in store.records if
        (label === nothing || record.label == label) &&
        _fingerprints_match(query, record.fingerprint)]
end

"""A conservative history-only suggestion; it never establishes mathematical eligibility."""
struct RouteAdvice
    matching_records::Int
    attempt_counts::Dict{Symbol, Int}
    success_counts::Dict{Symbol, Int}
    lower_confidence::Dict{Symbol, Float64}
    recommended_route::Union{Nothing, Symbol}
    selection_reason::Symbol
    preconditioner_reuse::Symbol
end

function _validate_history_policy(policy::HistoryPolicy)
    policy.min_samples > 0 || throw(ArgumentError("history min_samples must be positive"))
    policy.confidence_z >= 0 || throw(ArgumentError("history confidence_z must be nonnegative"))
    policy.exploration in (:off, :least_tried) ||
        throw(ArgumentError("history exploration must be :off or :least_tried"))
    return policy
end

function _wilson_lower_bound(successes::Int, attempts::Int, z::Float64)
    attempts == 0 && return 0.0
    proportion = successes / attempts
    z_squared = z^2
    denominator = 1 + z_squared / attempts
    center = proportion + z_squared / (2 * attempts)
    radius = z * sqrt(proportion * (1 - proportion) / attempts + z_squared / (4 * attempts^2))
    return max(0.0, (center - radius) / denominator)
end

function route_advice(store::HistoryStore, problem::AdaptiveLinearProblem,
        candidates::AbstractVector{Symbol}; profile::FingerprintProfile=FingerprintProfile(),
        conditioning_policy::ConditioningPolicy=ConditioningPolicy(),
        history_policy::HistoryPolicy=HistoryPolicy())
    _validate_history_policy(history_policy)
    query = fingerprint(problem; profile=profile, conditioning_policy=conditioning_policy)
    records = similar_records(store, query; label=problem.label)
    attempts = Dict{Symbol, Int}()
    successes = Dict{Symbol, Int}()
    for record in records
        record.route in candidates || continue
        attempts[record.route] = get(attempts, record.route, 0) + 1
        record.status in (Success, FallbackSuccess) &&
            (successes[record.route] = get(successes, record.route, 0) + 1)
    end
    confidence = Dict(route => _wilson_lower_bound(get(successes, route, 0),
        get(attempts, route, 0), history_policy.confidence_z) for route in candidates)
    recommended = nothing
    best_confidence = -1.0
    for route in candidates
        attempts_for_route = get(attempts, route, 0)
        if attempts_for_route >= history_policy.min_samples &&
           (recommended === nothing || confidence[route] > best_confidence)
            recommended = route
            best_confidence = confidence[route]
        end
    end
    selection_reason = recommended === nothing ? :insufficient_evidence : :confidence_ranked
    if recommended === nothing && history_policy.exploration == :least_tried && !isempty(candidates)
        recommended = first(candidates)
        for route in candidates
            get(attempts, route, 0) < get(attempts, recommended, 0) && (recommended = route)
        end
        selection_reason = :controlled_exploration
    end
    preconditioner = problem.preconditioner
    preconditioner_name = preconditioner === nothing ? :none : preconditioner.name
    reusable = preconditioner_name != :none && any(record -> record.status in (Success, FallbackSuccess) &&
        hasproperty(record, :preconditioner) && record.preconditioner == preconditioner_name, records)
    return RouteAdvice(length(records), attempts, successes, confidence, recommended, selection_reason,
        preconditioner_name == :none ? :not_applicable :
        (reusable ? :reuse_candidate : :no_reuse_evidence))
end

"""Persist only explicitly selected in-memory history; load files from trusted local paths only."""
function save_history(path::AbstractString, store::HistoryStore)
    open(path, "w") do io
        serialize(io, (format_version=1, records=store.records))
    end
    return path
end

function load_history!(store::HistoryStore, path::AbstractString)
    payload = open(deserialize, path)
    payload.format_version == 1 || throw(ArgumentError("unsupported history format"))
    for record in payload.records
        record isa NamedTuple || throw(ArgumentError("history contains an invalid record"))
        _record!(store, record)
    end
    return store
end

function _trace_payload(iteration, policy::TelemetryPolicy)
    policy.level in (:trace, :diagnostic) || return nothing
    iteration === nothing && return (residuals=Float64[], truncated=false)
    limit = policy.budget.max_trace_samples
    limit == 0 && return (residuals=Float64[], truncated=!isempty(iteration.residual_history))
    sampled = iteration.residual_history[1:policy.sample_every:end]
    truncated = length(sampled) > limit
    return (residuals=sampled[1:min(end, limit)], truncated=truncated)
end
