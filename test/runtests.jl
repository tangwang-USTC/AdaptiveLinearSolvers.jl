using Test
using AdaptiveLinearSolvers

@testset "AdaptiveLinearSolvers v0.1.2" begin
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

    generic = solve(A, b)
    @test generic.route == :lu
    @test isapprox(A * generic.x, b; rtol=1e-12)

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
