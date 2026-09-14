module AdaptiveLinearSolvers

using LinearAlgebra
using SparseArrays
using Krylov
using Serialization
using IterativeSolvers
import LinearSolve

# Optional external packages for preconditioner builders
const _HAS_INCOMPLETE_LU = try
    @eval using IncompleteLU
    true
catch
    @warn "IncompleteLU.jl not available — ILU preconditioner falls back to ILU(0)"
    false
end
const _HAS_ALGEBRAIC_MULTIGRID = try
    @eval using AlgebraicMultigrid
    true
catch
    @warn "AlgebraicMultigrid.jl not available — AMG preconditioner unavailable"
    false
end

include("types.jl")
include("operators.jl")
include("policy.jl")
include("iterative.jl")
include("preconditioners.jl")
include("krylov_backend.jl")
include("iterativesolvers_backend.jl")
include("linearsolve_backend.jl")
include("diagnostics.jl")
include("telemetry.jl")
include("backends.jl")
include("planning.jl")
include("routing.jl")

export AdaptiveLinearProblem, AdaptiveLinearSolution, Auto, Prefer, Lock, Forbid,
       RoutePolicy, MathematicalContract, PropertyEvidence, EvidenceLevel,
       Unknown, Suspected, Claimed, Certified, Proved, ConditioningInfo, SpectralInfo,
       DiagnosticBudget, TelemetryBudget, HistoryPolicy, ResourceBudget, BackendPolicy, BackendCapability,
       BackendSelection, ResourceAssessment, ConditioningPolicy, ConditioningAssessment, NumericalDiagnosis,
       OperatorApplicationBudgetExceeded, CountingOperator,
       PreconditionerContract, PreconditionerBuildPolicy, PreconditionerCache,
       PreconditionerBuildReport, build_preconditioner, with_preconditioner,
       clear_preconditioner_cache!, IterativeMethodCapability, IterativeQualification,
       iterative_capability, qualify_iterative,
       MatrixFreeOperator, BlockLayout, BlockOperator, blockrange,
       SolveStatus, Success, FallbackSuccess, QualificationRejected,
       NumericalFailure, BudgetTerminated, ResidualPolicy, RouteCertificate,
       IterationControl, IterationReport,
       RouteCapability, EligibilityDecision, LayerDecision, RoutePlan, plan,
       assess_conditioning, diagnose,
       TelemetryPolicy, OutputRequest, FingerprintProfile, MatrixFingerprint, RouteAdvice,
       HistoryStore, fingerprint, similar_records, route_advice, save_history, load_history!, solve,
       backend_capabilities, select_backend, assess_resources

end
