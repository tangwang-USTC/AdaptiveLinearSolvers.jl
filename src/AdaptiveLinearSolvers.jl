module AdaptiveLinearSolvers

using LinearAlgebra
using SparseArrays
using Krylov

include("types.jl")
include("operators.jl")
include("policy.jl")
include("iterative.jl")
include("krylov_backend.jl")
include("planning.jl")
include("telemetry.jl")
include("routing.jl")

export AdaptiveLinearProblem, AdaptiveLinearSolution, Auto, Prefer, Lock, Forbid,
       RoutePolicy, MathematicalContract, PropertyEvidence, EvidenceLevel,
       Unknown, Suspected, Claimed, Certified, Proved, ConditioningInfo,
       PreconditionerContract, IterativeMethodCapability, IterativeQualification,
       iterative_capability, qualify_iterative,
       MatrixFreeOperator, BlockLayout, BlockOperator, blockrange,
       SolveStatus, Success, FallbackSuccess, QualificationRejected,
       NumericalFailure, BudgetTerminated, ResidualPolicy, RouteCertificate,
       IterationControl, IterationReport,
       RouteCapability, EligibilityDecision, LayerDecision, RoutePlan, plan,
       TelemetryPolicy, OutputRequest, FingerprintProfile, MatrixFingerprint,
       HistoryStore, solve

end
