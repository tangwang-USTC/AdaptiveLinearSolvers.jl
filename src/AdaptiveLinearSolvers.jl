module AdaptiveLinearSolvers

using LinearAlgebra
using SparseArrays

include("types.jl")
include("policy.jl")
include("telemetry.jl")
include("routing.jl")

export AdaptiveLinearProblem, AdaptiveLinearSolution, Auto, Prefer, Lock, Forbid,
       RoutePolicy, MathematicalContract, PropertyEvidence, EvidenceLevel,
       Unknown, Suspected, Claimed, Certified, Proved, ConditioningInfo,
       TelemetryPolicy, OutputRequest, FingerprintProfile, MatrixFingerprint,
       HistoryStore, solve

end
