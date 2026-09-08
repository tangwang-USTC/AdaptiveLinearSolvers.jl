@testset "telemetry" begin
    A = standard_matrix()
    b = standard_rhs()
    contract = standard_spd_contract()

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
    @test AdaptiveLinearSolvers._validate_telemetry(TelemetryPolicy(level=:trace)).level == :trace
end
