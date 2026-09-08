@testset "operators" begin
    b = standard_rhs()
    matrix_free = MatrixFreeOperator(2, 2, (y, x) -> begin
        y[1] = 4 * x[1] + x[2]
        y[2] = x[1] + 3 * x[2]
        nothing
    end)
    matrix_free_problem = AdaptiveLinearProblem(matrix_free, b)
    matrix_free_plan = plan(matrix_free_problem)
    @test matrix_free_plan.execution_routes == [:gmres]
    @test all(decision -> !decision.eligible,
        matrix_free_plan.eligibility[1:4])
    matrix_free_result = solve(matrix_free_problem;
        iteration_control=IterationControl(max_iterations=20))
    @test matrix_free_result.status == Success
    @test matrix_free_result.route == :gmres

    locked_direct_matrix_free = plan(matrix_free_problem, RoutePolicy(direct=Lock(:lu)))
    @test isempty(locked_direct_matrix_free.planned_routes)

    layout = BlockLayout([1, 1])
    blocks = Matrix{Any}(undef, 2, 2)
    blocks[1, 1] = reshape([4.0], 1, 1)
    blocks[1, 2] = reshape([1.0], 1, 1)
    blocks[2, 1] = reshape([1.0], 1, 1)
    blocks[2, 2] = reshape([3.0], 1, 1)
    block_operator = BlockOperator(layout, blocks)
    block_problem = AdaptiveLinearProblem(block_operator, b)
    @test blockrange(layout, 2) == 2:2
    @test plan(block_problem).execution_routes == [:gmres]
    block_result = solve(block_problem;
        iteration_control=IterationControl(max_iterations=20))
    @test block_result.status == Success
    @test isapprox(block_operator * block_result.x, b; rtol=1e-12)
end
