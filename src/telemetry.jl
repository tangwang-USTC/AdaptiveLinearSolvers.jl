"""Request static fields from the current call; this does not create a fingerprint."""
Base.@kwdef struct OutputRequest
    route::Bool = false
    residual_ratio::Bool = false
    conditioning::Bool = false
end

"""Select labels needed for history-driven routing advice."""
Base.@kwdef struct FingerprintProfile
    representation::Bool = true
    size_band::Bool = true
    structure::Bool = true
    conditioning::Bool = true
    execution::Bool = true
end

"""Telemetry is off by default. Version 0.1.1 supports only `:off`, `:basic`, and `:fingerprint`."""
Base.@kwdef struct TelemetryPolicy
    level::Symbol = :off
    output::OutputRequest = OutputRequest()
    fingerprint::FingerprintProfile = FingerprintProfile()
    emit_on::Tuple{Vararg{Symbol}} = ()
end

function _validate_telemetry(policy::TelemetryPolicy)
    policy.level in (:off, :basic, :fingerprint) ||
        throw(ArgumentError("telemetry level $(policy.level) is planned after v0.1.1; use :off, :basic, or :fingerprint"))
    return policy
end

Base.@kwdef struct MatrixFingerprint
    representation::Symbol
    size_band::Symbol
    structure::Symbol
    conditioning::Symbol
    execution::Symbol = :serial_cpu
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

function _conditioning_tag(info::Union{Nothing, ConditioningInfo})
    info === nothing && return :unknown
    info.estimate === nothing && return :unknown
    info.estimate < 1e6 && return :moderate
    info.estimate < 1e12 && return :ill_conditioned
    return :severely_ill_conditioned
end

function _fingerprint(problem::AdaptiveLinearProblem)
    return MatrixFingerprint(
        representation=issparse(problem.A) ? :sparse_explicit : :dense_explicit,
        size_band=_size_band(problem.A),
        structure=_structure_tag(problem.contract),
        conditioning=_conditioning_tag(problem.conditioning),
    )
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
