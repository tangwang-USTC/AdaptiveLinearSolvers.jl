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

    qr_solution = solve(A, b; policy=RoutePolicy(direct=Lock(:qr)))
    @test qr_solution.route == :qr
    @test qr_solution.status == Success
    @test isapprox(A * qr_solution.x, b; rtol=1e-12)

    svd_solution = solve(A, b; policy=RoutePolicy(direct=Lock(:svd)))
    @test svd_solution.route == :svd
    @test svd_solution.status == Success
    @test isapprox(A * svd_solution.x, b; rtol=1e-12)

    direct_solution = solve(A, b; policy=RoutePolicy(direct=Lock(:direct)))
    @test direct_solution.route == :direct
    @test direct_solution.status == Success
    @test isapprox(A * direct_solution.x, b; rtol=1e-12)

    # 稀疏矩阵扩展：多种稀疏模式
    sparse_spd = sprand(10, 10, 0.3) + I
    sparse_spd = (sparse_spd + sparse_spd') / 2 + 10I  # 强制对角占优 SPD
    sparse_b = rand(10)
    sparse_contract = MathematicalContract(
        square=PropertyEvidence(Certified; source=:caller),
        hermitian=PropertyEvidence(Certified; source=:caller),
        positive_definite=PropertyEvidence(Certified; source=:caller),
    )
    sparse_spd_result = solve(AdaptiveLinearProblem(sparse_spd, sparse_b; contract=sparse_contract))
    @test sparse_spd_result.status == Success
    @test isapprox(sparse_spd * sparse_spd_result.x, sparse_b; rtol=1e-10)

    # 稀疏矩形矩阵 → QR（b 在列空间内可精确求解；确保每列非零）
    rect_x_true = rand(5)
    dense_rect = rand(8, 5)
    sparse_rect = sparse(dense_rect)
    rect_b_exact = sparse_rect * rect_x_true
    rect_result = solve(sparse_rect, rect_b_exact;
        policy=RoutePolicy(direct=Lock(:qr)),
        residual_policy=ResidualPolicy(absolute_tolerance=1e-10, relative_tolerance=1e-10))
    @test rect_result.status == Success
    @test rect_result.route == :qr
    @test isapprox(rect_x_true, rect_result.x; rtol=1e-8)

    # 对称不定矩阵 → Bunch-Kaufman（Hermitian 非 SPD）
    sym_indef = [1.0 2.0; 2.0 1.0]  # 特征值 3, -1
    sym_b = [1.0, 2.0]
    hermitian_contract = MathematicalContract(
        square=PropertyEvidence(Certified; source=:caller),
        hermitian=PropertyEvidence(Certified; source=:caller))
    bk_result = solve(AdaptiveLinearProblem(sym_indef, sym_b; contract=hermitian_contract);
        policy=RoutePolicy(direct=Lock(:bunchkaufman)))
    @test bk_result.status == Success
    @test bk_result.route == :bunchkaufman
    @test isapprox(sym_indef * bk_result.x, sym_b; rtol=1e-12)

    # SPD 矩阵上 Bunch-Kaufman 也应当可用
    spd_contract = standard_spd_contract()
    bk_spd = solve(AdaptiveLinearProblem(A, b; contract=spd_contract);
        policy=RoutePolicy(direct=Lock(:bunchkaufman)))
    @test bk_spd.status == Success

    # 非对称矩阵 → Bunch-Kaufman 拒绝（无 Hermitian 证据）
    bk_reject = solve([1.0 2.0; 3.0 4.0], b;
        policy=RoutePolicy(direct=Lock(:bunchkaufman)))
    @test bk_reject.status == QualificationRejected
end
