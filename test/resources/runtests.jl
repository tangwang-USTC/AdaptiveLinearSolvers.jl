@testset "resources" begin
    A = standard_matrix()
    b = standard_rhs()
    problem = AdaptiveLinearProblem(A, b)

    direct_backend = select_backend(problem, :lu)
    @test direct_backend.available
    @test direct_backend.name == :stdlib

    krylov_backend = select_backend(problem, :gmres)
    @test krylov_backend.available
    @test krylov_backend.name == :krylov

    iterative_solvers_backend = select_backend(problem, :gmres,
        BackendPolicy(requested=:iterativesolvers))
    @test iterative_solvers_backend.available
    @test iterative_solvers_backend.name == :iterativesolvers

    linearsolve_backend = select_backend(problem, :gmres,
        BackendPolicy(requested=:linearsolve))
    @test linearsolve_backend.available
    @test linearsolve_backend.name == :linearsolve

    unavailable_backend = select_backend(problem, :gmres,
        BackendPolicy(requested=:petsc))
    @test !unavailable_backend.available
    @test unavailable_backend.reason == :backend_adapter_not_implemented

    memory_limited = ResourceBudget(max_memory_bytes=1)
    assessment = assess_resources(problem, :lu, memory_limited)
    @test !assessment.eligible
    @test assessment.reason == :memory_budget_exceeded

    denied = solve(problem; resource_budget=memory_limited)
    @test denied.status == BudgetTerminated

    iteration_limited = solve(A, b;
        policy=RoutePolicy(iterative=Lock(:gmres)),
        iteration_control=IterationControl(max_iterations=20),
        resource_budget=ResourceBudget(max_iterations=1))
    @test iteration_limited.status == BudgetTerminated

    counted = solve(A, b;
        policy=RoutePolicy(iterative=Lock(:gmres)),
        iteration_control=IterationControl(max_iterations=20),
        resource_budget=ResourceBudget(max_operator_applications=4))
    @test counted.status == Success
    @test counted.iteration.operator_applications !== nothing
    @test counted.iteration.operator_applications <= 4

    alternative = solve(A, b;
        policy=RoutePolicy(iterative=Lock(:gmres)),
        backend_policy=BackendPolicy(requested=:iterativesolvers),
        iteration_control=IterationControl(max_iterations=20))
    @test alternative.status == Success
    @test alternative.route == :gmres

    ecosystem = solve(A, b;
        policy=RoutePolicy(iterative=Lock(:gmres)),
        backend_policy=BackendPolicy(requested=:linearsolve),
        iteration_control=IterationControl(max_iterations=20))
    @test ecosystem.status == Success
    @test ecosystem.route == :gmres

    # 边界值：未指定后端 → auto 自动选
    auto_backend = select_backend(problem, :cg)
    @test auto_backend.available
    @test auto_backend.name == :krylov

    # 边界值：未知后端名
    unknown_backend = select_backend(problem, :cg, BackendPolicy(requested=:nonexistent))
    @test !unknown_backend.available
    @test unknown_backend.reason == :unknown_backend

    # 边界值：未实现后端的 GPU 执行模式
    gpu_request = select_backend(problem, :cg,
        BackendPolicy(requested=:cuda))
    @test !gpu_request.available
    @test gpu_request.reason == :no_implemented_backend_for_route

    # 边界值：Full ResourceBudget 默认值
    default_budget = ResourceBudget()
    @test isinf(default_budget.max_seconds)
    @test default_budget.max_iterations == 0
    @test default_budget.max_memory_bytes == 0

    # 边界值：资源评估时未指定内存预算
    no_memory_check = assess_resources(problem, :lu, ResourceBudget())
    @test no_memory_check.eligible
    @test no_memory_check.reason == :unbounded_memory
end
