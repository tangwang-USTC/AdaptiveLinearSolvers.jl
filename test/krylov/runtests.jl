@testset "krylov" begin
    A = standard_matrix()
    b = standard_rhs()

    gmres_result = solve(A, b;
        policy=RoutePolicy(iterative=Lock(:gmres)),
        iteration_control=IterationControl(max_iterations=20, record_history=true))
    @test gmres_result.status == Success
    @test gmres_result.route == :gmres
    @test gmres_result.iteration.converged
    @test gmres_result.iteration.iterations > 0

    fgmres_result = solve(A, b;
        policy=RoutePolicy(iterative=Lock(:fgmres)),
        iteration_control=IterationControl(max_iterations=20))
    @test fgmres_result.status == Success
    @test fgmres_result.iteration.method == :fgmres
end
