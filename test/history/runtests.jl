@testset "history" begin
    A = standard_matrix()
    b = standard_rhs()
    profile = FingerprintProfile(size_band=false, structure=false,
        conditioning=false, execution=false)
    problem = AdaptiveLinearProblem(A, b)
    compact = fingerprint(problem; profile=profile)
    @test compact.representation == :dense_explicit
    @test compact.size_band === nothing
    @test compact.structure === nothing
    @test compact.conditioning === nothing
    @test compact.execution === nothing

    history = HistoryStore(4)
    solve(A, b; telemetry=TelemetryPolicy(level=:fingerprint, fingerprint=profile),
        history=history)
    advice = route_advice(history, problem, [:cholesky, :lu]; profile=profile)
    @test advice.matching_records == 1
    @test advice.recommended_route == :lu
    @test advice.preconditioner_reuse == :not_applicable
    @test first(plan(problem; history=history, fingerprint_profile=profile).execution_routes) == :lu

    traced = solve(A, b;
        policy=RoutePolicy(iterative=Lock(:gmres)),
        iteration_control=IterationControl(max_iterations=20, record_history=true),
        telemetry=TelemetryPolicy(level=:trace,
            budget=TelemetryBudget(max_trace_samples=2)))
    @test traced.status == Success
    @test length(traced.telemetry.trace.residuals) <= 2
end
