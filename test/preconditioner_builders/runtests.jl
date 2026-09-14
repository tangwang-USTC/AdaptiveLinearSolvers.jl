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

    # ─── IC(0) preconditioner (sparse SPD) ─────────────────────────────────────────

    @testset "IC(0) preconditioner" begin
        # 50×50 sparse SPD matrix — diagonal-dominant random SPD
        n_ic = 50
        A_spd = sprand(n_ic, n_ic, 0.15) + I
        A_spd = (A_spd + A_spd') / 2 + n_ic * I  # strongly diagonally dominant SPD
        b_ic = rand(n_ic)
        spd_contract_ic = MathematicalContract(
            square=PropertyEvidence(Certified; source=:caller),
            hermitian=PropertyEvidence(Certified; source=:caller),
            positive_definite=PropertyEvidence(Certified; source=:caller),
        )
        problem_ic = AdaptiveLinearProblem(A_spd, b_ic; contract=spd_contract_ic,
            label=:ic_test, matrix_version=:ic_v1)

        # Build
        ic_built = build_preconditioner(problem_ic; policy=PreconditionerBuildPolicy(kind=:ic))
        @test ic_built.reason == :built
        @test ic_built.contract.name == :ic
        @test ic_built.contract.operator isa AdaptiveLinearSolvers.ICZeroPreconditioner

        # Solve with CG
        ic_solved = solve(with_preconditioner(problem_ic, ic_built.contract);
            policy=RoutePolicy(iterative=Lock(:cg)),
            iteration_control=IterationControl(max_iterations=100))
        @test ic_solved.status == Success
        @test ic_solved.route == :cg
        @test isapprox(A_spd * ic_solved.x, b_ic; rtol=1e-6)
    end

    # IC(0) rejection: dense matrix
    @testset "IC(0) rejection: dense matrix" begin
        dense_contract = MathematicalContract(
            square=PropertyEvidence(Certified; source=:caller),
            hermitian=PropertyEvidence(Certified; source=:caller),
            positive_definite=PropertyEvidence(Certified; source=:caller),
        )
        dense_problem = AdaptiveLinearProblem(A, b; contract=dense_contract)
        rejected = build_preconditioner(dense_problem; policy=PreconditionerBuildPolicy(kind=:ic))
        @test rejected.reason == :sparse_matrix_required
        @test rejected.contract === nothing
    end

    # IC(0) rejection: no SPD contract
    @testset "IC(0) rejection: no SPD contract" begin
        n_ic = 20
        A_spd = sprand(n_ic, n_ic, 0.2) + I
        A_spd = (A_spd + A_spd') / 2 + n_ic * I
        no_contract_problem = AdaptiveLinearProblem(A_spd, rand(n_ic);
            label=:ic_nospd)
        rejected = build_preconditioner(no_contract_problem;
            policy=PreconditionerBuildPolicy(kind=:ic))
        @test rejected.reason == :spd_contract_required_for_ic
    end

    # ─── ILU preconditioner (sparse square) ─────────────────────────────────────────

    @testset "ILU preconditioner" begin
        n_ilu = 50
        # Random sparse non-symmetric matrix with diagonal dominance
        A_ilu = sprand(n_ilu, n_ilu, 0.15) + 2 * I
        b_ilu = rand(n_ilu)

        ilu_built = build_preconditioner(AdaptiveLinearProblem(A_ilu, b_ilu);
            policy=PreconditionerBuildPolicy(kind=:ilu))
        @test ilu_built.reason == :built
        @test ilu_built.contract.name == :ilu

        # Solve with GMRES (no Hermitian/SPD required)
        ilu_solved = solve(with_preconditioner(AdaptiveLinearProblem(A_ilu, b_ilu), ilu_built.contract);
            policy=RoutePolicy(iterative=Lock(:gmres)),
            iteration_control=IterationControl(max_iterations=100))
        @test ilu_solved.status == Success
        @test ilu_solved.route == :gmres
        @test isapprox(A_ilu * ilu_solved.x, b_ilu; rtol=1e-6)
    end

    # ILU rejection: dense matrix
    @testset "ILU rejection: dense matrix" begin
        rejected = build_preconditioner(AdaptiveLinearProblem(A, b);
            policy=PreconditionerBuildPolicy(kind=:ilu))
        @test rejected.reason == :sparse_matrix_required
    end

    # ILU rejection: non-square
    @testset "ILU rejection: non-square" begin
        rect_sparse = sparse(rand(8, 5))
        rejected = build_preconditioner(AdaptiveLinearProblem(rect_sparse, rand(8));
            policy=PreconditionerBuildPolicy(kind=:ilu))
        @test rejected.reason == :nonsquare_operator
    end

    # ─── AMG preconditioner (sparse SPD) ────────────────────────────────────────────

    @testset "AMG preconditioner" begin
        n_amg = 50
        A_amg = sprand(n_amg, n_amg, 0.15) + I
        A_amg = (A_amg + A_amg') / 2 + n_amg * I  # SPD
        b_amg = rand(n_amg)
        spd_amg = MathematicalContract(
            square=PropertyEvidence(Certified; source=:caller),
            hermitian=PropertyEvidence(Certified; source=:caller),
            positive_definite=PropertyEvidence(Certified; source=:caller),
        )

        amg_built = build_preconditioner(
            AdaptiveLinearProblem(A_amg, b_amg; contract=spd_amg);
            policy=PreconditionerBuildPolicy(kind=:amg))
        @test amg_built.reason == :built
        @test amg_built.contract.name == :amg

        amg_solved = solve(
            with_preconditioner(AdaptiveLinearProblem(A_amg, b_amg; contract=spd_amg), amg_built.contract);
            policy=RoutePolicy(iterative=Lock(:cg)),
            iteration_control=IterationControl(max_iterations=100))
        @test amg_solved.status == Success
        @test amg_solved.route == :cg
        @test isapprox(A_amg * amg_solved.x, b_amg; rtol=1e-6)

        # Ruge-Stuben variant
        amg_rs_built = build_preconditioner(
            AdaptiveLinearProblem(A_amg, b_amg; contract=spd_amg);
            policy=PreconditionerBuildPolicy(kind=:amg, amg_type=:ruge_stuben))
        @test amg_rs_built.reason == :built
        @test amg_rs_built.contract.name == :amg
    end

    # AMG rejection: dense matrix
    @testset "AMG rejection: dense matrix" begin
        rejected = build_preconditioner(AdaptiveLinearProblem(A, b);
            policy=PreconditionerBuildPolicy(kind=:amg))
        @test rejected.reason == :sparse_matrix_required
    end

    # ─── Policy validation ──────────────────────────────────────────────────────────

    @testset "policy validation" begin
        @test_throws ArgumentError AdaptiveLinearSolvers._validate_preconditioner_policy(
            PreconditionerBuildPolicy(kind=:nonexistent))
        @test_throws ArgumentError AdaptiveLinearSolvers._validate_preconditioner_policy(
            PreconditionerBuildPolicy(diagonal_tolerance=-1.0))
        @test_throws ArgumentError AdaptiveLinearSolvers._validate_preconditioner_policy(
            PreconditionerBuildPolicy(ilu_droptol=-0.1))
        @test_throws ArgumentError AdaptiveLinearSolvers._validate_preconditioner_policy(
            PreconditionerBuildPolicy(ilu_tau=-0.1))
        @test_throws ArgumentError AdaptiveLinearSolvers._validate_preconditioner_policy(
            PreconditionerBuildPolicy(amg_type=:invalid))
        @test AdaptiveLinearSolvers._validate_preconditioner_policy(
            PreconditionerBuildPolicy(kind=:ilu)).kind == :ilu
        @test AdaptiveLinearSolvers._validate_preconditioner_policy(
            PreconditionerBuildPolicy(kind=:ic)).kind == :ic
        @test AdaptiveLinearSolvers._validate_preconditioner_policy(
            PreconditionerBuildPolicy(kind=:amg)).kind == :amg
    end

    # ─── ICZeroPreconditioner ldiv! directly ────────────────────────────────────────

    @testset "ICZeroPreconditioner forward/backward substitution" begin
        # Small SPD system where we can verify the solve
        n_small = 10
        A_small = sparse(Diagonal(rand(1.0:2.0, n_small)) + sprand(n_small, n_small, 0.3))
        A_small = (A_small + A_small') / 2 + n_small * I  # ensure SPD
        L = AdaptiveLinearSolvers._ic_zero_factor(A_small)
        M = AdaptiveLinearSolvers.ICZeroPreconditioner(L)

        # Verify M \ x is correct
        x = rand(n_small)
        y = similar(x)
        AdaptiveLinearSolvers.ldiv!(y, M, x)
        # M * y should ≈ x (i.e., L * L' * y ≈ x)
        @test isapprox(M.L * (M.L' * y), x; rtol=1e-10)
    end
end