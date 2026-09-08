@enum EvidenceLevel begin
    Unknown
    Suspected
    Claimed
    Certified
    Proved
end

"""Evidence for one mathematical property used by route qualification."""
struct PropertyEvidence
    level::EvidenceLevel
    source::Symbol
end

PropertyEvidence(level::EvidenceLevel=Unknown; source::Symbol=:unspecified) =
    PropertyEvidence(level, source)

"""Mathematical facts supplied by the caller; unknown is intentionally distinct from false."""
Base.@kwdef struct MathematicalContract
    square::PropertyEvidence = PropertyEvidence()
    hermitian::PropertyEvidence = PropertyEvidence()
    positive_definite::PropertyEvidence = PropertyEvidence()
    nonsingular::PropertyEvidence = PropertyEvidence()
    rank_deficient::PropertyEvidence = PropertyEvidence()
end

"""Optional conditioning evidence obtained outside the router."""
Base.@kwdef struct ConditioningInfo
    estimate::Union{Nothing, Float64} = nothing
    source::Symbol = :unspecified
    reliable::Bool = false
end

"""A linear system plus optional mathematical and operational evidence."""
Base.@kwdef struct AdaptiveLinearProblem{TA, TB, TC, TI}
    A::TA
    b::TB
    contract::TC = MathematicalContract()
    conditioning::TI = nothing
    label::Symbol = :anonymous
end

function AdaptiveLinearProblem(A, b;
        contract::MathematicalContract=MathematicalContract(),
        conditioning::Union{Nothing, ConditioningInfo}=nothing,
        label::Symbol=:anonymous)
    return AdaptiveLinearProblem(A, b, contract, conditioning, label)
end

"""Result of one solve, including only telemetry enabled by the selected policy."""
struct AdaptiveLinearSolution{TX, TT, TH}
    x::TX
    status::Symbol
    route::Symbol
    residual_ratio::Float64
    telemetry::TT
    history::TH
end

_certified_or_proved(evidence::PropertyEvidence) =
    evidence.level == Certified || evidence.level == Proved
