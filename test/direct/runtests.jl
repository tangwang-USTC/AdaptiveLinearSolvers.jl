@testset "direct" begin
    A = standard_matrix()
    b = standard_rhs()
    contract = standard_spd_contract()

    result = solve(AdaptiveLinearProblem(A, b; contract=contract))
    @test result.route == :cholesky
    @test result.status == Success
    @test result.residual_ratio < 1e-12

    generic = solve(A, b)
    @test generic.route == :lu
    @test isapprox(A * generic.x, b; rtol=1e-12)

    sparse_result = solve(sparse(A), b)
    @test sparse_result.status == Success

    zero_rhs = solve(A, zeros(2); residual_policy=ResidualPolicy(absolute_tolerance=1e-12))
    @test zero_rhs.status == Success
    @test zero_rhs.residual_ratio === nothing
end
