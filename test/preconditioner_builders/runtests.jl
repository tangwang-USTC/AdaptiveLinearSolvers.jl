@testset "preconditioner builders" begin
    A = standard_matrix()
    b = standard_rhs()
    contract = standard_spd_contract()
    problem = AdaptiveLinearProblem(A, b; contract=contract, label=:vfp_step,
        matrix_version=:assembly_1)
    cache = PreconditionerCache(2)

    built = build_preconditioner(problem; cache=cache)
    @test built.reason == :built
    @test !built.cache_hit
    @test built.contract.name == :jacobi
    @test qualify_iterative(with_preconditioner(problem, built.contract), :cg).eligible

    reused = build_preconditioner(problem; cache=cache)
    @test reused.cache_hit
    @test reused.contract === built.contract

    changed_version = AdaptiveLinearProblem(A, b; contract=contract, label=:vfp_step,
        matrix_version=:assembly_2)
    rebuilt = build_preconditioner(changed_version; cache=cache)
    @test !rebuilt.cache_hit

    solved = solve(with_preconditioner(problem, built.contract);
        policy=RoutePolicy(iterative=Lock(:cg)),
        iteration_control=IterationControl(max_iterations=20))
    @test solved.status == Success
    @test solved.route == :cg

    identity_built = build_preconditioner(problem; policy=PreconditionerBuildPolicy(kind=:identity))
    @test identity_built.reason == :built
    @test identity_built.contract.name == :identity
    @test identity_built.contract.operator ≈ Diagonal(ones(2))

    identity_solved = solve(with_preconditioner(problem, identity_built.contract);
        policy=RoutePolicy(iterative=Lock(:cg)),
        iteration_control=IterationControl(max_iterations=20))
    @test identity_solved.status == Success
end
