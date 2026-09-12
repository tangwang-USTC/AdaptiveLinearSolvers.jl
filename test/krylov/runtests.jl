@testset "krylov" begin
    A = standard_matrix()
    b = standard_rhs()
    spd = standard_spd_contract()
    hermitian_contract = MathematicalContract(
        square=PropertyEvidence(Certified; source=:caller),
        hermitian=PropertyEvidence(Certified; source=:caller),
    )

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

    cg_result = solve(AdaptiveLinearProblem(A, b; contract=spd);
        policy=RoutePolicy(iterative=Lock(:cg)),
        iteration_control=IterationControl(max_iterations=20, record_history=true))
    @test cg_result.status == Success
    @test cg_result.route == :cg
    @test cg_result.iteration.converged
    @test isapprox(A * cg_result.x, b; rtol=1e-12)

    minres_result = solve(AdaptiveLinearProblem(A, b; contract=hermitian_contract);
        policy=RoutePolicy(iterative=Lock(:minres)),
        iteration_control=IterationControl(max_iterations=20, record_history=true))
    @test minres_result.status == Success
    @test minres_result.route == :minres
    @test minres_result.iteration.converged
    @test isapprox(A * minres_result.x, b; rtol=1e-6)

    bicgstab_result = solve(A, b;
        policy=RoutePolicy(iterative=Lock(:bicgstab)),
        iteration_control=IterationControl(max_iterations=40, record_history=true))
    @test bicgstab_result.status == Success
    @test bicgstab_result.route == :bicgstab
    @test bicgstab_result.iteration.converged
    @test isapprox(A * bicgstab_result.x, b; rtol=1e-6)
end

@testset "iteration control validation" begin
    @test_throws ArgumentError AdaptiveLinearSolvers._validate_iteration_control(IterationControl(max_iterations=-1))
    @test_throws ArgumentError AdaptiveLinearSolvers._validate_iteration_control(IterationControl(max_seconds=0))
    @test AdaptiveLinearSolvers._validate_iteration_control(IterationControl(max_iterations=0)).max_iterations == 0
    @test isinf(AdaptiveLinearSolvers._validate_iteration_control(IterationControl()).max_seconds)
end
