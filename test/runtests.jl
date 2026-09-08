using Test
using SparseArrays
using AdaptiveLinearSolvers

@testset "AdaptiveLinearSolvers v0.0.2" begin
    A = [4.0 1.0; 1.0 3.0]
    b = [1.0, 2.0]
    contract = MathematicalContract(
        square=PropertyEvidence(Certified; source=:caller),
        hermitian=PropertyEvidence(Certified; source=:caller),
        positive_definite=PropertyEvidence(Certified; source=:caller),
    )
    result = solve(AdaptiveLinearProblem(A, b; contract=contract))
    @test result.route == :cholesky
    @test result.status == Success
    @test result.residual_ratio < 1e-12

    direct_plan = plan(AdaptiveLinearProblem(A, b; contract=contract),
        RoutePolicy(family=Lock(:direct), direct=Lock(:cholesky)))
    @test direct_plan.execution_routes == [:cholesky]
    @test all(decision -> decision.accepted, direct_plan.layer_decisions[1:2])

    unqualified_plan = plan(AdaptiveLinearProblem(A, b), RoutePolicy(direct=Lock(:cholesky)))
    @test isempty(unqualified_plan.execution_routes)
    @test only(unqualified_plan.eligibility).reason == :hermitian_evidence_insufficient

    iterative_plan = plan(AdaptiveLinearProblem(A, b), RoutePolicy(iterative=Lock(:gmres)))
    @test iterative_plan.execution_routes == [:gmres]
    @test iterative_plan.planned_routes == [:gmres]
    @test isempty(iterative_plan.unavailable_routes)

    gmres_result = solve(A, b;
        policy=RoutePolicy(iterative=Lock(:gmres)),
        iteration_control=IterationControl(max_iterations=20, record_history=true))
    @test gmres_result.status == Success
    @test gmres_result.route == :gmres
    @test gmres_result.iteration.converged
    @test gmres_result.iteration.iterations > 0

    generic = solve(A, b)
    @test generic.route == :lu
    @test isapprox(A * generic.x, b; rtol=1e-12)

    sparse_result = solve(sparse(A), b)
    @test sparse_result.status == Success

    zero_rhs = solve(A, zeros(2); residual_policy=ResidualPolicy(absolute_tolerance=1e-12))
    @test zero_rhs.status == Success
    @test zero_rhs.residual_ratio === nothing

    rank_deficient = [1.0 0.0; 0.0 0.0]
    rank_contract = MathematicalContract(rank_deficient=PropertyEvidence(Certified; source=:caller))
    rank_plan = plan(AdaptiveLinearProblem(rank_deficient, [1.0, 0.0]; contract=rank_contract))
    @test first(rank_plan.execution_routes) == :svd

    fixed_preconditioner = PreconditionerContract(
        name=:ilu,
        fixed_within_solve=PropertyEvidence(Certified; source=:caller),
        linear_within_solve=PropertyEvidence(Certified; source=:caller),
    )
    variable_preconditioner = PreconditionerContract(name=:adaptive_ilu)
    fixed_problem = AdaptiveLinearProblem(A, b; contract=contract, preconditioner=fixed_preconditioner)
    variable_problem = AdaptiveLinearProblem(A, b; contract=contract, preconditioner=variable_preconditioner)
    @test qualify_iterative(fixed_problem, :cg).eligible
    @test qualify_iterative(fixed_problem, :gmres).eligible
    @test qualify_iterative(variable_problem, :gmres).reason == :variable_or_unknown_preconditioner_requires_fgmres
    @test qualify_iterative(variable_problem, :fgmres).eligible
    @test !qualify_iterative(AdaptiveLinearProblem(A, b), :minres).eligible

    flexible_plan = plan(variable_problem, RoutePolicy(iterative=Prefer(:gmres)))
    @test flexible_plan.planned_routes[1] == :fgmres
    @test flexible_plan.unavailable_routes == [:fgmres]

    fgmres_result = solve(A, b;
        policy=RoutePolicy(iterative=Lock(:fgmres)),
        iteration_control=IterationControl(max_iterations=20))
    @test fgmres_result.status == Success
    @test fgmres_result.iteration.method == :fgmres

    history = HistoryStore(1)
    telemetry = TelemetryPolicy(level=:fingerprint,
        output=OutputRequest(route=true, residual_ratio=true), emit_on=(:success,))
    traced = solve(A, b; telemetry=telemetry, history=history)
    @test traced.telemetry.fingerprint.representation == :dense_explicit
    @test length(history.records) == 1

    audited = solve(AdaptiveLinearProblem(A, b; contract=contract),
        telemetry=TelemetryPolicy(level=:basic, output=OutputRequest(certificate=true)))
    @test audited.certificate.selected_route == :cholesky
    @test audited.certificate.residual_accepted

    rejected = solve(A, b; policy=RoutePolicy(direct=Lock(:cholesky)))
    @test rejected.status == QualificationRejected
    @test_throws ArgumentError AdaptiveLinearSolvers._validate_telemetry(TelemetryPolicy(level=:trace))
end
