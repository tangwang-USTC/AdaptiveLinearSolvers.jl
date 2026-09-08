@testset "planning" begin
    A = standard_matrix()
    b = standard_rhs()
    contract = standard_spd_contract()

    direct_plan = plan(AdaptiveLinearProblem(A, b; contract=contract),
        RoutePolicy(family=Lock(:direct), direct=Lock(:cholesky)))
    @test direct_plan.execution_routes == [:cholesky]
    @test all(decision -> decision.accepted, direct_plan.layer_decisions[1:2])

    unqualified_plan = plan(AdaptiveLinearProblem(A, b), RoutePolicy(direct=Lock(:cholesky)))
    @test isempty(unqualified_plan.execution_routes)
    @test only(unqualified_plan.eligibility).reason == :hermitian_evidence_insufficient

    iterative_plan = plan(AdaptiveLinearProblem(A, b), RoutePolicy(iterative=Lock(:gmres)))
    @test iterative_plan.execution_routes == [:gmres]
    @test iterative_plan.planned_routes == [:gmres]
    @test isempty(iterative_plan.unavailable_routes)

    rank_deficient = [1.0 0.0; 0.0 0.0]
    rank_contract = MathematicalContract(rank_deficient=PropertyEvidence(Certified; source=:caller))
    rank_plan = plan(AdaptiveLinearProblem(rank_deficient, [1.0, 0.0]; contract=rank_contract))
    @test first(rank_plan.execution_routes) == :svd
end
