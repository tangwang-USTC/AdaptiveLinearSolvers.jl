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
end
