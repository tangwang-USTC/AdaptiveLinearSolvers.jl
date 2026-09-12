@testset "fallback" begin
    A = standard_matrix()
    b = standard_rhs()

    # 回退到直接法：迭代 Lock(:bicgstab) 无 contract，bicgstab 收敛失败 → 回退到直接法
    fallback_result = solve(A, b;
        policy=RoutePolicy(iterative=Prefer(:bicgstab)),
        iteration_control=IterationControl(max_iterations=40))
    @test fallback_result.status == FallbackSuccess || fallback_result.status == Success
    @test isapprox(A * fallback_result.x, b; rtol=1e-12)

    # 全路线被拒绝：Lock(:cholesky) 但无 Hermitian 证据
    rejected = solve(A, b; policy=RoutePolicy(direct=Lock(:cholesky)))
    @test rejected.status == QualificationRejected
    @test rejected.route === nothing
    @test rejected.certificate === nothing

    # 内存预算导致所有路线被终止
    tiny_memory = solve(A, b; resource_budget=ResourceBudget(max_memory_bytes=1))
    @test tiny_memory.status == BudgetTerminated
    @test tiny_memory.x === nothing || isempty(tiny_memory.x)

    # 迭代路线锁定时预算终止（telemetry 默认 :off → certificate=nothing，仅验证 status）
    budget_exceeded = solve(A, b;
        policy=RoutePolicy(iterative=Lock(:gmres)),
        iteration_control=IterationControl(max_iterations=20),
        resource_budget=ResourceBudget(max_iterations=1))
    @test budget_exceeded.status == BudgetTerminated

    # 无 contract 时优先走 LU，SPD contract 时走 Cholesky
    spd = standard_spd_contract()
    spd_result = solve(AdaptiveLinearProblem(A, b; contract=spd))
    @test spd_result.route == :cholesky
    @test spd_result.status == Success
end