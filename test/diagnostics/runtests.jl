@testset "diagnostics" begin
    A = standard_matrix()
    b = standard_rhs()
    versioned_information = ConditioningInfo(
        estimate=1e14,
        metric=:kappa_2,
        operator=:original,
        matrix_version=:assembly_17,
        source=:caller_estimate,
        reliable=true,
        evidence=:external,
    )
    versioned_problem = AdaptiveLinearProblem(A, b;
        conditioning=versioned_information, matrix_version=:assembly_17)
    accepted = assess_conditioning(versioned_problem)
    @test accepted.accepted
    @test accepted.state == :near_rank_deficient
    @test first(plan(versioned_problem).execution_routes) == :svd

    stale_problem = AdaptiveLinearProblem(A, b;
        conditioning=versioned_information, matrix_version=:assembly_18)
    stale = assess_conditioning(stale_problem)
    @test !stale.accepted
    @test stale.reason == :matrix_version_mismatch
    @test stale.state == :stale

    no_estimate = diagnose(AdaptiveLinearProblem(A, b; matrix_version=:assembly_17))
    @test !no_estimate.estimate_performed
    @test no_estimate.conditioning.reason == :no_conditioning_information

    estimate_policy = ConditioningPolicy(estimation=:full,
        budget=DiagnosticBudget(max_seconds=1.0, max_matrix_dimension=4))
    estimated = diagnose(AdaptiveLinearProblem(A, b; matrix_version=:assembly_17);
        policy=estimate_policy)
    @test estimated.estimate_performed
    @test estimated.conditioning.information.metric == :kappa_2
    @test estimated.conditioning.information.matrix_version == :assembly_17

    variable_preconditioner = PreconditionerContract(name=:adaptive_ilu,
        operator=Diagonal(ones(2)))
    stagnated_report = IterationReport(:fgmres, 10, false, "maximum iterations reached",
        0.01, [1.0, 0.98, 0.97])
    stagnated = diagnose(AdaptiveLinearProblem(A, b;
        preconditioner=variable_preconditioner);
        iteration=stagnated_report, status=NumericalFailure)
    @test stagnated.iteration_state == :stagnated
    @test stagnated.preconditioner_state == :suspected_failure
    @test stagnated.overall_state == :iterative_stagnation

    breakdown_report = IterationReport(:gmres, 3, false, "Krylov breakdown", 0.01, Float64[])
    breakdown = diagnose(AdaptiveLinearProblem(A, b); iteration=breakdown_report,
        status=NumericalFailure)
    @test breakdown.iteration_state == :breakdown
    @test breakdown.overall_state == :iterative_breakdown
end
