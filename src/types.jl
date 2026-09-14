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

"""Conditioning information with explicit metric, operator, provenance, and matrix version."""
Base.@kwdef struct ConditioningInfo
    estimate::Union{Nothing, Float64} = nothing
    metric::Symbol = :unknown
    operator::Symbol = :original
    matrix_version::Any = nothing
    source::Symbol = :unspecified
    reliable::Bool = false
    evidence::Symbol = :unspecified
end

"""Approximate spectral interval obtained under an explicit Hermitian-operator contract."""
Base.@kwdef struct SpectralInfo
    lambda_min::Union{Nothing, Float64} = nothing
    lambda_max::Union{Nothing, Float64} = nothing
    method::Symbol = :unknown
    operator::Symbol = :original
    matrix_version::Any = nothing
    source::Symbol = :unspecified
    steps::Int = 0
    reliable::Bool = false
end

"""Hard limits for opt-in numerical diagnostics; zero dimension means no estimate is permitted."""
Base.@kwdef struct DiagnosticBudget
    max_seconds::Float64 = 0.0
    max_operator_applications::Int = 0
    max_matrix_dimension::Int = 0
end

"""Caps for opt-in telemetry sampling; zero trace samples disables trace retention."""
Base.@kwdef struct TelemetryBudget
    max_trace_samples::Int = 0
    max_extra_operator_applications::Int = 0
    max_seconds::Float64 = 0.0
end

"""Statistical safeguards for history-based ranking; exploration is explicit and deterministic."""
Base.@kwdef struct HistoryPolicy
    min_samples::Int = 3
    confidence_z::Float64 = 1.96
    exploration::Symbol = :off
end

"""Execution limits shared by route selection and backend invocation; zero integer limits mean no extra cap."""
Base.@kwdef struct ResourceBudget
    max_seconds::Float64 = Inf
    max_iterations::Int = 0
    max_memory_bytes::Int = 0
    max_operator_applications::Int = 0
    execution::Symbol = :serial_cpu
end

"""Backend selection is explicit; `:auto` chooses only registered, implemented capabilities."""
Base.@kwdef struct BackendPolicy
    requested::Symbol = :auto
    allow_fallback::Bool = true
end

"""Declared backend coverage, independent of whether an optional package is installed."""
struct BackendCapability
    name::Symbol
    package::Symbol
    routes::Vector{Symbol}
    execution_modes::Vector{Symbol}
    implemented::Bool
end

"""Result of matching a route to a backend without invoking numerical work."""
struct BackendSelection
    name::Union{Nothing, Symbol}
    available::Bool
    reason::Symbol
end

"""Conservative resource admission result; unknown memory is never silently treated as bounded."""
struct ResourceAssessment
    eligible::Bool
    reason::Symbol
    estimated_memory_bytes::Union{Nothing, Int}
end

"""Raised internally when a counted iterative operator reaches its explicit call limit."""
struct OperatorApplicationBudgetExceeded <: Exception
    limit::Int
end

Base.showerror(io::IO, error::OperatorApplicationBudgetExceeded) =
    print(io, "operator application budget exceeded (limit=$(error.limit))")

"""Policy for validating caller data and optionally estimating conditioning after route planning."""
Base.@kwdef struct ConditioningPolicy
    use_external::Bool = true
    require_matching_version::Bool = true
    require_reliable_external::Bool = true
    estimation::Symbol = :none
    spectral_estimation::Symbol = :none
    spectral_steps::Int = 12
    budget::DiagnosticBudget = DiagnosticBudget()
end

"""Validation result for supplied or explicitly estimated conditioning information."""
struct ConditioningAssessment
    information::Union{Nothing, ConditioningInfo}
    accepted::Bool
    reason::Symbol
    state::Symbol
end

"""Separated numerical diagnosis; no state is silently promoted to a mathematical qualification."""
struct NumericalDiagnosis
    conditioning::ConditioningAssessment
    spectral::Union{Nothing, SpectralInfo}
    spectral_state::Symbol
    iteration_state::Symbol
    preconditioner_state::Symbol
    overall_state::Symbol
    estimate_performed::Bool
    elapsed_seconds::Float64
    notes::Vector{String}
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

"""Construction policy for preconditioners implemented by this package."""
Base.@kwdef struct PreconditionerBuildPolicy
    kind::Symbol = :jacobi
    diagonal_tolerance::Float64 = sqrt(eps(Float64))
    reuse::Bool = true
    ilu_droptol::Float64 = 0.01
    ilu_tau::Float64 = 0.001
    amg_type::Symbol = :smoothed_aggregation
end

"""Bounded in-memory cache keyed by an explicit matrix version and construction policy."""
mutable struct PreconditionerCache
    capacity::Int
    entries::Dict{Any, PreconditionerContract}
    order::Vector{Any}
end

PreconditionerCache(capacity::Integer=16) =
    PreconditionerCache(Int(capacity), Dict{Any, PreconditionerContract}(), Any[])

"""Result of an explicit preconditioner build or cache reuse request."""
struct PreconditionerBuildReport
    contract::Union{Nothing, PreconditionerContract}
    cache_key::Any
    cache_hit::Bool
    build_seconds::Float64
    reason::Symbol
end

"""A linear system plus optional mathematical and operational evidence."""
Base.@kwdef struct AdaptiveLinearProblem{TA, TB, TC, TI, TP, TV}
    A::TA
    b::TB
    contract::TC = MathematicalContract()
    conditioning::TI = nothing
    preconditioner::TP = nothing
    label::Symbol = :anonymous
    matrix_version::TV = nothing
end

function AdaptiveLinearProblem(A, b;
        contract::MathematicalContract=MathematicalContract(),
        conditioning::Union{Nothing, ConditioningInfo}=nothing,
        preconditioner::Union{Nothing, PreconditionerContract}=nothing,
        label::Symbol=:anonymous, matrix_version=nothing)
    return AdaptiveLinearProblem(A, b, contract, conditioning, preconditioner, label,
        matrix_version)
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
