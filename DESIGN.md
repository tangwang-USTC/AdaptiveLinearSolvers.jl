# AdaptiveLinearSolvers.jl 软件设计

## 术语与缩略词

| 缩略词 | 首次定义 |
|---|---|
| API | Application Programming Interface，应用程序编程接口 |
| LU | Lower-Upper factorization，下三角-上三角分解 |
| QR | Orthogonal-Upper Triangular factorization，正交-上三角分解 |
| SVD | Singular Value Decomposition，奇异值分解 |
| CG | Conjugate Gradient，共轭梯度法 |
| MINRES | Minimum Residual，最小残量法 |
| QLP | QLP decomposition，QLP 分解 |
| GMRES | Generalized Minimum Residual，广义最小残量法 |
| FGMRES | Flexible Generalized Minimum Residual，柔性广义最小残量法 |
| LSQR | Least-Squares QR，最小二乘 QR 法 |
| LSMR | Least-Squares Minimum Residual，最小二乘最小残量法 |
| AMG | Algebraic Multigrid，代数多重网格 |
| ILU | Incomplete LU，不完全 LU 分解 |
| PETSc | Portable, Extensible Toolkit for Scientific Computation，可移植可扩展科学计算工具包 |
| MPI | Message Passing Interface，消息传递接口 |
| GPU | Graphics Processing Unit，图形处理器 |

## 1. 设计原则

1. 数学核心不依赖物理领域语义。
2. 路由资格先于性能启发式。
3. 用户可以锁定任意决策层，其余层保持自动。
4. 后端是可选能力，不能污染核心 API。
5. 每次结果必须带可机器读取和可人工审查的路线证书。

## 2. 公共对象

```julia
struct AdaptiveLinearProblem{A,B,X,C}
    A::A
    b::B
    x0::X
    contract::C
end

struct AdaptiveLinearSolution{X,C,H}
    x::X
    certificate::C
    history::H
end
```

`AdaptiveLinearProblem` 不直接包装某个后端的问题对象。它保留原始线性算子、右端项、初值、数学契约和缓存键；执行时再映射到 Julia 标准库、Krylov 或 `LinearSolve.jl`。

## 3. 分层路线覆盖

`ForcedRoute` 不表示必须锁定完整路线。它是分层覆盖策略，每一层都可取 `Auto()`、`Prefer(value)`、`Lock(value)` 或 `Forbid(value)`。

```julia
policy = RoutePolicy(
    family          = Lock(:krylov),
    method          = Lock(:fgmres),
    preconditioner  = Auto(),
    backend         = Prefer(:krylovjl),
    precision       = Auto(),
    tolerances      = Auto(),
    fallback        = Auto(),
)

solution = solve(problem, policy)
```

可覆盖的层级为：

| 层级 | 示例 | 未锁定时的行为 |
|---|---|---|
| 问题表述 | 显式矩阵、矩阵自由、最小二乘 | 从 `A` 与维度推断 |
| 求解器族 | 直接法、Krylov、最小二乘 | 按资格与预算选择 |
| 具体方法 | CG、MINRES、FGMRES、QR | 在选定族内选择 |
| 预条件器 | 无、Jacobi、ILU、AMG、用户算子 | 按结构、稀疏度和后端能力选择 |
| 后端 | 标准库、Krylov、LinearSolve、PETSc | 选择已加载且兼容的后端 |
| 精度 | 原精度、混合精度、提升精度 | 依据数值风险和硬件能力选择 |
| 容差和预算 | 相对/绝对容差、最大迭代、内存 | 根据问题尺度与默认策略选择 |
| 回退 | 允许、禁止、限定集合 | 自动生成有限回退图 |

优先级固定为：不可绕过的数学资格 > 显式 `Lock` > `Prefer` > 自动策略。`Lock(:cg)` 若没有 CG 所需证据，默认返回 `IneligibleRoute`；只有显式 `unsafe=true` 才允许执行，并在证书中记录该越权。

## 4. 核心模块

```text
src/
  AdaptiveLinearSolvers.jl
  problem.jl             问题对象与公共 solve API
  contracts.jl           数学契约、证据与零空间
  inspection.jl          显式矩阵结构检查
  eligibility.jl         算法资格规则
  policy.jl              Auto/Prefer/Lock/Forbid 与合并
  planning.jl            路线图、预算和回退图
  cache.jl               分解、预条件器和模式缓存
  verification.jl        残差、后向误差与状态
  certificate.jl         可序列化的路线证书
  backends/stdlib.jl     LinearAlgebra 与 SparseArrays
ext/
  AdaptiveLinearSolversKrylovExt.jl
  AdaptiveLinearSolversLinearSolveExt.jl
  AdaptiveLinearSolversAMGExt.jl
  AdaptiveLinearSolversILUExt.jl
  AdaptiveLinearSolversPETScExt.jl
  AdaptiveLinearSolversCUDAExt.jl
```

## 5. 后端能力模型

每个扩展注册纯描述性的能力：支持的标量类型、设备、显式/矩阵自由表示、所需结构、预条件器语义和缓存种类。计划器只读取能力，不直接依赖外部包对象。

| 扩展 | 初始职责 |
|---|---|
| 标准库后端 | LU、Cholesky、Bunch-Kaufman、QR、SVD、稀疏直接法 |
| Krylov 扩展 | CG、MINRES-QLP、GMRES、FGMRES、BiCGStab、LSQR/LSMR |
| LinearSolve 扩展 | `LinearProblem` 映射、`init`/`solve!` 缓存和后端适配 |
| AMG 扩展 | 已认证 SPD 系统的 AMG 预条件器 |
| ILU 扩展 | 一般稀疏系统的 ILU 预条件器 |
| PETSc 扩展 | MPI、字段分裂、Schur 补和大型分布式问题 |

Krylov 路线以 `Krylov.jl` 为主，因为其同时覆盖 Hermitian、非 Hermitian、最小二乘和矩阵自由问题[[1](#ref-1)]。`LinearSolve.jl` 的价值在于统一后端与缓存，而其默认选路不替代本项目的资格层[[2](#ref-2)]。

## 6. 路由与执行

```text
inspect(problem)
  -> merge(contract, observed_evidence)
  -> apply(policy_overrides)
  -> reject(ineligible_routes)
  -> rank(eligible_routes)
  -> build(primary + bounded_fallbacks)
  -> execute
  -> verify
  -> emit(certificate)
```

路由器不应以完整条件数作为大规模问题的常规前置计算。对迭代法，执行监控应使用残差下降率、停滞、breakdown、内存、时间预算和预条件器构造成本。

预条件器只有在同一次 `solve(A, b)` 内保持固定线性算子时，才可选择 GMRES；若其在迭代中变化、是非线性的，或该性质未知，规划器必须选择 FGMRES。能力模型通过 `fixed_within_solve` 与 `linear_within_solve` 表达该条件；数学理由、适用边界和例外见[理论基础的固定与可变预条件器章节](THEORY.md#fixed-and-variable-preconditioners)。

## 7. 缓存

缓存分为三类：

- 模式缓存：稀疏结构分析、块分区和重排序。
- 数值缓存：LU/Cholesky/QR 分解及 AMG 层次。
- 迭代缓存：初值、Krylov 工作区和预条件器状态。

调用方必须提供可复用键或数值版本号。核心不通过扫描完整矩阵内容猜测“矩阵是否改变”。

## 8. 验证、状态与测试

标准状态为 `Converged`、`Stagnated`、`Breakdown`、`IneligibleRoute`、`NumericallySuspect` 和 `BudgetExceeded`。`Converged` 必须由统一残差验收器确认，不能仅转发外部后端的成功标志。

第一批测试只包含数学实例：稠密/稀疏 SPD、Hermitian 不定、非对称非正规、矩形、秩亏、矩阵自由、固定与可变预条件器，以及人工制造的算法资格违例。

## 参考文献列表

[1]. <a id="ref-1"></a> JuliaSmoothOptimizers. 2026. [*Matrix-free operators*](https://github.com/JuliaSmoothOptimizers/Krylov.jl/blob/main/docs/src/matrix_free.md). Krylov.jl Documentation.<br>

[2]. <a id="ref-2"></a> SciML. 2026. [*Linear System Solvers*](https://docs.sciml.ai/LinearSolve/stable/solvers/solvers/). LinearSolve.jl Documentation.<br>
