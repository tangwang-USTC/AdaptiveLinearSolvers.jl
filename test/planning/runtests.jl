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

    # 组合策略：direct Lock(:qr)、无 SPD 证据 → QR 唯一可执行
    qr_plan = plan(AdaptiveLinearProblem(A, b), RoutePolicy(direct=Lock(:qr)))
    @test qr_plan.execution_routes == [:qr]
    @test length(qr_plan.planned_routes) == 1

    # 组合策略：iterative=Lock(:cg) + SPD 无证据 → 空执行路线
    cg_no_spd = plan(AdaptiveLinearProblem(A, b), RoutePolicy(iterative=Lock(:cg)))
    @test isempty(cg_no_spd.execution_routes)
    @test cg_no_spd.eligibility[1].reason == :hermitian_evidence_insufficient

    # 组合策略：Prefer(:gmres) + 可变预条件器 → FGMRES 替代 GMRES
    variable_pc = PreconditionerContract(name=:adaptive, operator=Diagonal(ones(2)))
    prefer_gmres = plan(AdaptiveLinearProblem(A, b; preconditioner=variable_pc),
        RoutePolicy(iterative=Prefer(:gmres)))
    @test prefer_gmres.planned_routes[1] == :fgmres
    @test :lu in prefer_gmres.planned_routes

    # 组合策略：Lock(:iterative) family + Lock(:cg) 无 SPD → 空
    family_cg = plan(AdaptiveLinearProblem(A, b),
        RoutePolicy(family=Lock(:iterative), iterative=Lock(:cg)))
    @test isempty(family_cg.execution_routes)
end
