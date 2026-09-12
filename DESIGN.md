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
| rcond | Reciprocal Condition Number，倒数条件数 |

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

条件信息是数学契约的可选字段，而不是每次调用都触发的前置计算：

```julia
struct ConditioningInfo
    estimate::Union{Nothing,Float64}
    metric::Symbol          # :kappa_1, :kappa_2, :kappa_inf, :rcond
    operator::Symbol        # :original, :left_preconditioned, :right_preconditioned
    matrix_version
    source::Symbol
    reliable::Bool
    evidence::Symbol        # :exact, :estimated, :external, :qualitative
end
```

`ConditioningPolicy` 控制是否使用外部信息、是否要求版本/可信度匹配、是否允许 `:cheap` 或 `:full` 估计，以及诊断预算。`plan` 只验证外部信息，保持无数值执行副作用；`diagnose` 才在调用方明确授权预算时估计。版本失配、不可靠或无效信息不会改变路线排序。`NumericalDiagnosis` 将条件状态、迭代状态、预条件器状态与综合状态分开返回；其按需输出由 `OutputRequest(diagnostics=true)` 或路线证书承载。完整语义、估计方法与未知条件数时的路线见[理论基础的规模、资源与条件信息章节](THEORY.md#scale-and-conditioning)。

```julia
struct TelemetryPolicy{S,K}
    level::Symbol               # :off, :basic, :fingerprint, :trace, :diagnostic
    output_request
    fingerprint_profile
    emit_on
    sample_every::Int
    capture_true_residual::Bool
    history_store::S
    family_key::K
end
```

默认 `TelemetryPolicy` 为 `level=:off, emit_on=()`。`AdaptiveLinearSolution` 仍返回解和轻量 `SolveStatus`，但 `RouteCertificate`、残差轨迹和 `HistoryStore` 写入均不生成。单次调用需要固定字段时，使用 `level=:basic` 与 `output_request`，不构造指纹。希望积累历史建议或按问题标签选择输出字段时，显式设置 `level=:fingerprint, emit_on=(:new_fingerprint, :failure, :fallback, :periodic)` 与 `fingerprint_profile`；`trace` 与 `diagnostic` 必须由调用方显式请求，并受时间、内存和额外算子应用预算约束。完整的输出语义、标签指纹和历史安全边界见[理论基础的可观测性、历史与性能预算章节](THEORY.md#51-可观测性历史与性能预算)。

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

路由器不应以完整条件数作为大规模问题的常规前置计算，但必须接受版本匹配的外部条件信息，并在诊断预算允许时执行按需估计。对迭代法，执行监控应使用残差下降率、停滞、breakdown、内存、时间预算和预条件器构造成本；具体路线规则见[理论基础的规模、资源与条件信息章节](THEORY.md#scale-and-conditioning)。

`ResourceBudget` 统一收紧时间、迭代、内存和算子应用预算；`BackendPolicy` 只可选择已登记且已实现的后端。当前标准库直接法与 `Krylov.jl` 仅支持串行 CPU；GPU、MPI、PETSc、LinearSolve.jl 与 IterativeSolvers.jl 作为未实现适配器显式登记，不能被自动选择。后端覆盖、资源准入和后续接入契约见[计算资源与后端](docs/BACKENDS.md)。

`PreconditionerBuildPolicy`、`PreconditionerCache` 与 `build_preconditioner` 形成显式构造和版本键复用层；当前只实现 Jacobi 与恒等预条件器，并由 `with_preconditioner` 接入既有 Krylov 执行路径。缓存没有 `matrix_version` 时自动禁用，不得由相似指纹或历史成功替代数值版本匹配。构造器范围与 ILU、AMG、Schur、物理预条件器的后续接入边界见[预条件器构造与版本缓存](docs/PRECONDITIONERS.md)。

`ConditioningPolicy(spectral_estimation=:lanczos)` 在时间与算子调用预算明确时，可对已认证 Hermitian 的显式或矩阵自由算子执行 Lanczos 谱区间诊断。该结果只进入 `NumericalDiagnosis`，不能反向建立 Hermitian/正定资格或替代外部条件信息。完整假设与解释边界见[大规模条件与谱诊断](docs/SPECTRAL_DIAGNOSTICS.md)。

`HistoryStore` 只接收 `fingerprint` 及以上级别、且同时满足 `emit_on` 事件条件的 `SolveRecord`，并按 `family_key` 与标签指纹检索相似历史。它只能向计划器提供候选路线的排序分数、预条件器复用提示和诊断建议；资格门、用户 `Lock`/`Forbid` 和当前预算仍具有更高优先级。初始实现使用容量受限的内存存储，避免在默认求解路径中写入大型原始数据。

## 9. `0.0.3` 实施边界

首个 Julia 实现仅覆盖显式稠密和稀疏矩阵的直接路线：通用反斜杠、`LU`（Lower-Upper）分解、Cholesky 分解、`QR`（Orthogonal-Triangular）分解与 `SVD`（Singular Value Decomposition，奇异值分解）。自动路线只有在调用方以 `Certified` 或 `Proved` 给出 Hermitian 正定证据时才选 Cholesky；它不会从矩阵元素抽样推断该资格。

`RoutePolicy` 的 `family`、`direct`、`iterative`、`preconditioner` 与 `fallback` 均接受 `Auto`、`Prefer`、`Lock` 或 `Forbid`，因此可以只锁定一个层级，其余层级继续保持自动。`0.0.3` 会将未注册后端的迭代路线标记为不可执行，而不会将其交给直接执行器。

实现支持 `off`、`basic`、`fingerprint`、`trace` 与 `diagnostic` 遥测。`FingerprintProfile` 实际裁剪返回和匹配字段；`HistoryStore` 提供相似检索、成功率建议与显式的受信任本机持久化。历史建议只重排已合格的自动候选路线，不能覆盖用户策略或数学资格。`TelemetryBudget` 限制 trace 样本；条件估计仍须由独立的 `ConditioningPolicy` 显式授权。`SolveStatus` 区分成功、回退成功、资格拒绝、数值失败与预算终止；`RouteCertificate` 在请求时保存候选路线、资格证据、尝试历史、回退原因、残差验收和数值诊断。实现细节与持久化安全边界见[可观测性与历史自适应](docs/OBSERVABILITY.md)。

## 10. `0.0.3` 规划批次

当前开发批次将规划与执行分离。`plan(problem, policy)` 返回 `RoutePlan`，其中包含候选路线、逐路线 `EligibilityDecision`、各策略层的 `LayerDecision` 和可执行路线顺序；该函数不得分解矩阵、申请求解工作区或计算残差。`solve` 只执行 `RoutePlan.execution_routes`，并将计划判定写入按需 `RouteCertificate`。

在当前基线中，直接方法族及其 `direct` 层已可执行；`iterative`、`preconditioner` 与 `fallback` 层已进入独立判定模型，但尚无可执行能力。对这些层的 `Lock` 请求必须给出空计划和明确拒绝，对 `Prefer` 请求仅记录“当前不可用”并允许其他自动路线继续。只有该规划批次的测试、文档和验收条件整体完成后，才递增版本号。

## 11. `0.0.3` 迭代资格批次

本批次先定义迭代法与预条件器的统一数学资格接口，不调用迭代内核。`PreconditionerContract` 记录预条件器是否在一次 `solve` 内固定，以及是否在线性意义下固定；`qualify_iterative(problem, method)` 返回可审计的 `IterativeQualification`。

首批资格规则覆盖 CG（Conjugate Gradient，共轭梯度法）、MINRES（Minimum Residual，最小残量法）、GMRES（Generalized Minimal Residual，广义最小残量法）、FGMRES（Flexible Generalized Minimal Residual，柔性广义最小残量法）和 BiCGStab（Biconjugate Gradient Stabilized，稳定化双共轭梯度法）。CG 要求已认证 Hermitian 正定性，MINRES 要求已认证 Hermitian 性；对要求固定线性预条件器的方法，若调用方提供的预条件器在一次求解内可变、非线性或性质未知，资格判定必须拒绝并给出 `:variable_or_unknown_preconditioner_requires_fgmres`。FGMRES 为该情形的可行候选，但本批次不执行任何 Krylov 迭代。

## 12. `0.0.3` 迭代路线规划批次

`RoutePlan` 将路线分为 `candidate_routes`、`planned_routes`、`execution_routes` 与 `unavailable_routes`。`planned_routes` 只表示数学资格通过；`execution_routes` 还要求当前存在可调用的数值后端；`unavailable_routes` 使调用方和审计记录能够区分“数学上正确但尚未实现”与“数学资格不满足”。

当 `iterative=Lock(:gmres)` 时，规划器只保留 GMRES；若其预条件器语义不满足固定线性条件，则没有合格路线。当 `iterative=Prefer(:gmres)` 且预条件器可变、非线性或未知时，规划器保留 GMRES 的拒绝记录并将 FGMRES 放入合格计划。当前所有迭代路线仍属于 `unavailable_routes`；下一开发任务是接入首个可执行的迭代后端及其停止、失败和残差监控接口。

## 13. `0.0.3` Krylov 执行批次

本批次以 `Krylov.jl` 作为首个执行后端，接入 CG、MINRES、GMRES、FGMRES 与 BiCGStab。`IterationControl` 统一提供最大迭代数、墙钟时间、重启选项和残差历史记录开关；`IterationReport` 统一返回收敛标志、完成迭代数、后端状态文本、实际耗时和残差历史。`solve` 在后端未收敛且到达调用方迭代或时间上限时返回 `BudgetTerminated`，其余未收敛情形返回 `NumericalFailure`。

该后端批次最初只执行无预条件器迭代路线；可应用预条件器接口已在同一发布批次中补齐。

## 14. `0.0.3` 预条件器应用批次

`PreconditionerContract.operator` 表示预条件器逆作用 $P^{-1}$，而不是矩阵 $P$ 本身。调用方负责提供可被 Krylov 后端应用的线性算子或矩阵；本项目在此批次不构造 Jacobi、ILU（Incomplete LU，不完全 LU 分解）、代数多重网格或物理预条件器。

对固定线性预条件器，CG、MINRES、GMRES 与 BiCGStab 将该逆作用传入 Krylov 的左预条件参数 `M`。对 FGMRES，可变或非线性预条件器作为柔性右预条件参数 `N` 传入。CG 与 MINRES 在存在预条件器时还要求调用方以 `Certified` 或 `Proved` 证明该预条件器 Hermitian 正定；证据不足必须拒绝，而非凭名称或矩阵抽样猜测。没有 `operator` 的非空预条件器契约仍使路线保持不可执行。

## 15. `0.0.3` 矩阵自由算子批次

`MatrixFreeOperator(rows, cols, apply!)` 只要求调用方实现原位线性作用 $y=Ax$，并提供 `size`、`mul!` 与向量乘法语义，因此无需显式组装稠密或稀疏矩阵即可交给 Krylov 后端。矩阵自由的数学定义、显式稀疏矩阵的区别、Vlasov-Fokker-Planck/Maxwell 例子、残差语义和伴随作用限制见[理论文档 3.1 节](THEORY.md#matrix-free-operators)。

路由器对矩阵自由表示的 LU、Cholesky、QR 和 SVD 返回 `:matrix_free_direct_route_unavailable`；若调用方未锁定直接路线，且没有其他直接路线合格，规划器自动生成满足数学资格的迭代候选。当前接口尚未定义伴随作用 $A^\dagger x$、GPU、MPI 或分布式向量。

## 16. `0.0.3` 块耦合算子批次

`BlockLayout` 将扁平 Krylov 向量划分为连续字段切片，`BlockOperator(layout, blocks)` 按块累加显式或矩阵自由子块作用。它既可表达 Vlasov-Fokker-Planck 高阶矩与 Maxwell 场的全耦合线性化，也不会要求后端改变向量存储格式。数学定义、字段耦合例子及其与 Schur 补/字段分裂的关系见[理论文档 3.2 节](THEORY.md#block-coupled-operators)。

块算子当前仍作为矩阵自由整体路由到 Krylov；没有实现块三角、Schur 补、字段分裂预条件器或块伴随作用。后续任何物理预条件器必须显式声明所依赖的字段布局、近似和数学资格。

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
