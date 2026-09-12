@testset "spectral" begin
    A = standard_matrix()
    b = standard_rhs()
    contract = standard_spd_contract()
    problem = AdaptiveLinearProblem(A, b; contract=contract, matrix_version=:assembly_1)
    policy = ConditioningPolicy(spectral_estimation=:lanczos, spectral_steps=2,
        budget=DiagnosticBudget(max_seconds=1.0, max_operator_applications=2))
    diagnosis = diagnose(problem; policy=policy)
    @test diagnosis.spectral.method == :lanczos
    @test diagnosis.spectral.steps == 2
    @test diagnosis.spectral.lambda_min > 0
    @test diagnosis.spectral.lambda_max > diagnosis.spectral.lambda_min
    @test diagnosis.spectral_state == :estimated

    matrix_free = MatrixFreeOperator(2, 2, (y, x) -> begin
        y[1] = 4 * x[1] + x[2]
        y[2] = x[1] + 3 * x[2]
        nothing
    end)
    matrix_free_problem = AdaptiveLinearProblem(matrix_free, b; contract=contract)
    matrix_free_diagnosis = diagnose(matrix_free_problem; policy=policy)
    @test matrix_free_diagnosis.spectral.method == :lanczos

    unqualified = diagnose(AdaptiveLinearProblem(A, b); policy=policy)
    @test unqualified.spectral === nothing
    @test unqualified.spectral_state == :unavailable
end
