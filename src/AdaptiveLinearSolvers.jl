module AdaptiveLinearSolvers

using LinearAlgebra
using SparseArrays

include("types.jl")
include("policy.jl")
include("planning.jl")
include("telemetry.jl")
include("routing.jl")

export AdaptiveLinearProblem, AdaptiveLinearSolution, Auto, Prefer, Lock, Forbid,
       RoutePolicy, MathematicalContract, PropertyEvidence, EvidenceLevel,
       Unknown, Suspected, Claimed, Certified, Proved, ConditioningInfo,
       SolveStatus, Success, FallbackSuccess, QualificationRejected,
       NumericalFailure, BudgetTerminated, ResidualPolicy, RouteCertificate,
       RouteCapability, EligibilityDecision, LayerDecision, RoutePlan, plan,
       TelemetryPolicy, OutputRequest, FingerprintProfile, MatrixFingerprint,
       HistoryStore, solve

end
