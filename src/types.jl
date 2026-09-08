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

"""Mathematical semantics of a preconditioner supplied by the caller."""
Base.@kwdef struct PreconditionerContract
    name::Symbol = :none
    operator::Any = nothing
    fixed_within_solve::PropertyEvidence = PropertyEvidence()
    linear_within_solve::PropertyEvidence = PropertyEvidence()
    hermitian::PropertyEvidence = PropertyEvidence()
    positive_definite::PropertyEvidence = PropertyEvidence()
end

"""A linear system plus optional mathematical and operational evidence."""
Base.@kwdef struct AdaptiveLinearProblem{TA, TB, TC, TI, TP}
    A::TA
    b::TB
    contract::TC = MathematicalContract()
    conditioning::TI = nothing
    preconditioner::TP = nothing
    label::Symbol = :anonymous
end

function AdaptiveLinearProblem(A, b;
        contract::MathematicalContract=MathematicalContract(),
        conditioning::Union{Nothing, ConditioningInfo}=nothing,
        preconditioner::Union{Nothing, PreconditionerContract}=nothing,
        label::Symbol=:anonymous)
    return AdaptiveLinearProblem(A, b, contract, conditioning, preconditioner, label)
end

@enum SolveStatus begin
    Success
    FallbackSuccess
    QualificationRejected
    NumericalFailure
    BudgetTerminated
end

"""Residual acceptance policy; a zero right-hand side is assessed by absolute residual only."""
Base.@kwdef struct ResidualPolicy
    absolute_tolerance::Float64 = 0.0
    relative_tolerance::Float64 = sqrt(eps(Float64))
end

"""Backend-independent iteration and wall-clock limits for one Krylov solve."""
Base.@kwdef struct IterationControl
    max_iterations::Int = 0
    max_seconds::Float64 = Inf
    restart::Bool = false
    record_history::Bool = false
end

"""Result of one solve, including only telemetry enabled by the selected policy."""
struct AdaptiveLinearSolution{TX, TC, TI, TT, TH}
    x::TX
    status::SolveStatus
    route::Union{Nothing, Symbol}
    residual_ratio::Union{Nothing, Float64}
    certificate::TC
    iteration::TI
    telemetry::TT
    history::TH
end

_certified_or_proved(evidence::PropertyEvidence) =
    evidence.level == Certified || evidence.level == Proved

_status_event(status::SolveStatus) = status == Success ? :success :
    status == FallbackSuccess ? :fallback :
    status == QualificationRejected ? :qualification_rejected :
    status == NumericalFailure ? :numerical_failure : :budget_terminated
