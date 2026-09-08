using Test
using LinearAlgebra
using SparseArrays
using AdaptiveLinearSolvers

include("test_support.jl")

const DEFAULT_TEST_GROUPS = (
    "direct", "planning", "krylov", "preconditioner", "operators", "diagnostics", "history", "resources", "telemetry",
)
const VALID_TEST_TIERS = ("default", "nightly", "manual", "all")
const TEST_TIER = get(ENV, "ALS_TEST_TIER", "default")

TEST_TIER in VALID_TEST_TIERS || error(
    "Unknown ALS_TEST_TIER=$(repr(TEST_TIER)); choose default, nightly, manual, or all.")

include_test_group(group::AbstractString) = include(joinpath(@__DIR__, group, "runtests.jl"))

function include_optional_tier(tier::AbstractString)
    directory = joinpath(@__DIR__, tier)
    isdir(directory) || return
    for filename in sort(readdir(directory))
        endswith(filename, ".jl") && include(joinpath(directory, filename))
    end
end

@testset "AdaptiveLinearSolvers v0.0.4" begin
    foreach(include_test_group, DEFAULT_TEST_GROUPS)
    TEST_TIER in ("nightly", "all") && include_optional_tier("nightly")
    TEST_TIER in ("manual", "all") && include_optional_tier("manual")
end
