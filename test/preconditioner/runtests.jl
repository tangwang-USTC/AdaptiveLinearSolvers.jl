@testset "preconditioner" begin
    A = standard_matrix()
    b = standard_rhs()
    contract = standard_spd_contract()

    fixed_preconditioner = PreconditionerContract(
        name=:ilu,
        operator=Diagonal([0.25, 1 / 3]),
        fixed_within_solve=PropertyEvidence(Certified; source=:caller),
        linear_within_solve=PropertyEvidence(Certified; source=:caller),
        hermitian=PropertyEvidence(Certified; source=:caller),
        positive_definite=PropertyEvidence(Certified; source=:caller),
    )
    variable_preconditioner = PreconditionerContract(
        name=:adaptive_ilu,
        operator=Diagonal(ones(2)),
    )
    fixed_problem = AdaptiveLinearProblem(A, b; contract=contract, preconditioner=fixed_preconditioner)
    variable_problem = AdaptiveLinearProblem(A, b; contract=contract, preconditioner=variable_preconditioner)
    @test qualify_iterative(fixed_problem, :cg).eligible
    @test qualify_iterative(fixed_problem, :gmres).eligible
    @test qualify_iterative(variable_problem, :gmres).reason == :variable_or_unknown_preconditioner_requires_fgmres
    @test qualify_iterative(variable_problem, :fgmres).eligible
    @test !qualify_iterative(AdaptiveLinearProblem(A, b), :minres).eligible

    insufficient_preconditioner = PreconditionerContract(
        name=:unspecified_spd,
        operator=Diagonal(ones(2)),
        fixed_within_solve=PropertyEvidence(Certified; source=:caller),
        linear_within_solve=PropertyEvidence(Certified; source=:caller),
    )
    @test qualify_iterative(AdaptiveLinearProblem(A, b; contract=contract,
        preconditioner=insufficient_preconditioner), :cg).reason ==
        :preconditioner_hermitian_positive_definite_evidence_insufficient

    flexible_plan = plan(variable_problem, RoutePolicy(iterative=Prefer(:gmres)))
    @test flexible_plan.planned_routes[1] == :fgmres
    @test flexible_plan.execution_routes[1] == :fgmres
    @test isempty(flexible_plan.unavailable_routes)

    flexible_result = solve(variable_problem;
        policy=RoutePolicy(iterative=Prefer(:gmres)),
        iteration_control=IterationControl(max_iterations=20))
    @test flexible_result.status == Success
    @test flexible_result.route == :fgmres
end
