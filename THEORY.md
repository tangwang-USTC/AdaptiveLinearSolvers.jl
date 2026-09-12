# AdaptiveLinearSolvers.jl 理论基础

## 术语与缩略词

| 缩略词 | 首次定义 |
|---|---|
| LU | Lower-Upper factorization，下三角-上三角分解 |
| QR | Orthogonal-Upper Triangular factorization，正交-上三角分解 |
| SVD | Singular Value Decomposition，奇异值分解 |
| SPD | Symmetric Positive Definite，对称正定 |
| CG | Conjugate Gradient，共轭梯度法 |
| LDL | Lower-Diagonal-Lower-Transpose factorization，下三角-对角-下三角转置分解 |
| MINRES | Minimum Residual，最小残量法 |
| QLP | QLP decomposition，QLP 分解 |
| GMRES | Generalized Minimum Residual，广义最小残量法 |
| FGMRES | Flexible Generalized Minimum Residual，柔性广义最小残量法 |
| LSQR | Least-Squares QR，最小二乘 QR 法 |
| LSMR | Least-Squares Minimum Residual，最小二乘最小残量法 |
| AMG | Algebraic Multigrid，代数多重网格 |
| PETSc | Portable, Extensible Toolkit for Scientific Computation，可移植可扩展科学计算工具包 |
| GPU | Graphics Processing Unit，图形处理器 |
| LAPACK | Linear Algebra PACKage，线性代数软件包 |
| rcond | Reciprocal Condition Number，倒数条件数 |

## 1. 范围与问题定义

本项目处理线性问题

$$
A x = b,
$$

以及最小二乘、最小范数和带显式零空间约束的推广。核心层只接受数学对象：显式矩阵、块矩阵或实现了 `mul!` 的线性算子；它不依赖具体物理模型。

调用方可提供已证明的数学契约，例如 Hermitian 性、正定性、零空间、块分区和可复用性。求解器也可从显式矩阵得到有限的结构证据。所有路线选择必须区分“已证明/已认证”和“仅怀疑”。

Julia 标准库提供 LU、Cholesky、Bunch-Kaufman、QR、SVD、稀疏分解与可复用的 `factorize` 原语，构成直接法基线[[1](#ref-1)]。矩阵自由路线以 Krylov 方法为主，并要求算子实现必要的乘法与伴随操作[[2](#ref-2)]。

## 2. 数学契约与证据

每个性质都由值、证据等级和来源组成：

| 证据等级 | 含义 | 自动路由权限 |
|---|---|---|
| `proved` | 来自离散推导、定理或显式同构变换 | 可授权严格方法 |
| `certified` | 由当前问题的完整结构检查或数值分解成功确认 | 可授权严格方法 |
| `claimed` | 调用方声明，尚未由本次问题核验 | 仅 `GuidedRoute` 可使用 |
| `suspected` | 抽样、历史经验或近似判断 | 只能影响排序 |
| `unknown` | 没有信息 | 不授权依赖该性质的方法 |

例如，CG 的自动资格为：方阵、Hermitian、正定，且预条件化时预条件器满足相应正定要求。MINRES 需要 Hermitian 性；Bunch-Kaufman/LDL 需要对称或 Hermitian 结构。非对称或证据不足时，通用 Krylov 主路线是 GMRES；预条件器随迭代变化时必须使用 FGMRES，而非普通 GMRES。

“数值上接近对称”不是 Hermitian 的证明。对显式矩阵，`Symmetric` 或 `Hermitian` 包装、全矩阵检查与成功的 Cholesky 分解可形成可记录的证据；对矩阵自由算子，必须由调用方提供能被验证的伴随与内积契约。

## 3. 线性问题分类

| 问题类别 | 数学条件 | 初始路线 | 保守回退 |
|---|---|---|---|
| 小型稠密方阵 | 一般可逆 | 主元 LU | QR，SVD |
| SPD | 已认证 Hermitian 正定 | Cholesky 或 PCG | AMG-PCG，直接分解 |
| Hermitian 不定或奇异 | 已认证 Hermitian | LDL/Bunch-Kaufman，MINRES 或 MINRES-QLP | 直接分解或最小范数路线 |
| 一般非对称方阵 | 无额外结构 | 重启 GMRES | FGMRES，BiCGStab |
| 超定/欠定 | 矩形 | QR | SVD，LSQR/LSMR |
| 近秩亏 | 秩揭示证据 | 枢轴 QR | 截断 SVD |
| 矩阵自由 | `mul!` 可用 | Krylov | 依据伴随是否可用选择最小二乘路线 |

<a id="matrix-free-operators"></a>
### 3.1 矩阵自由算子

设 $X$ 和 $Y$ 为有限维向量空间，线性映射为

$$
\mathcal{A}: X \rightarrow Y, \qquad y = \mathcal{A}x.
$$

矩阵自由算子（matrix-free operator）指调用方能够计算向量作用 $y=\mathcal{A}x$，但不在内存中显式存储其相对于某组基的矩阵 $A$。这里“矩阵自由”不表示矩阵在数学上不存在；它表示路由器和求解器只能访问算子作用、维度及调用方提供的结构契约，而不能逐项读取 $A_{ij}$、进行稀疏模式分析或构造直接分解。Krylov 方法只反复需要该向量作用，因而适合此表示[[2](#ref-2)]。

显式稀疏矩阵与矩阵自由算子不同：前者仍保存非零位置和数值，可使用稀疏 LU、fill-in 分析或 ILU（Incomplete LU，不完全 LU 分解）；后者通常以离散算子组合、局部核、快速变换或雅可比向量积实现，主要成本是一次 $x\mapsto\mathcal{A}x$ 的计算。

当前接口要求原位作用：

```julia
A = MatrixFreeOperator(n, n, apply!)

function apply!(y, x)
    # overwrite y with A * x
    return y
end
```

`size(A)=(m,n)` 给出定义域和值域维度，`mul!(y, A, x)` 实现 $y\leftarrow\mathcal{A}x$。残差仍可按

$$
r=b-\mathcal{A}x
$$

计算，因此相对残差 $\lVert r\rVert/\lVert b\rVert$ 的定义不因不显式组装矩阵而改变；当 $b=0$ 时仍只使用绝对残差。

**Vlasov-Fokker-Planck 例子。** 在给定时间步或 Newton 线性化点 $u_0$，离散残量 $F(u)$ 的雅可比作用可直接实现为

$$
\mathcal{A}v \approx \frac{F(u_0+\epsilon v)-F(u_0)}{\epsilon},
$$

或由流动、碰撞和场耦合离散块逐项作用得到。这样避免形成高维相空间中的完整雅可比矩阵；代价是每次 Krylov 迭代需要一次残量或分块算子求值。差分近似是否足够线性、$\epsilon$ 如何选择以及边界条件如何处理，必须由物理离散层记录为契约，不能由线性代数路由器猜测。

**Maxwell 例子。** 对频域或隐式时间步的电场未知量，可将 curl-curl 型作用写为

$$
\mathcal{A}E = \nabla\times\left(\mu^{-1}\nabla\times E\right) + \sigma E,
$$

其中离散旋度、材料系数、导电项和边界条件在 `apply!` 内部依次应用。若该离散在指定内积下 Hermitian 或正定，调用方必须以证据形式明确提供；仅从连续方程形式或局部抽样不能授权 CG 或 MINRES。

矩阵自由表示的路线边界如下：

- LU、Cholesky、QR 和 SVD 需要显式元素或专门的隐式分解算法，当前直接后端不具备该资格。
- CG、MINRES、GMRES、FGMRES 和 BiCGStab 只需正向作用，但各自仍受 Hermitian 性、正定性和预条件器语义资格门约束。
- 最小二乘、条件数一范数估计和涉及 $A^\ast$ 的路线还需要伴随作用 $z\mapsto\mathcal{A}^\ast z$；当前 `MatrixFreeOperator` 未声明此能力，因此不得自动选择这些路线。

<a id="block-coupled-operators"></a>
### 3.2 块耦合算子

对多物理或多字段未知量，将全局向量按字段分块为

$$
x=
\begin{bmatrix}
x_1\\
\vdots\\
x_p
\end{bmatrix},
\qquad
\mathcal{A}=
\begin{bmatrix}
\mathcal{A}_{11} & \cdots & \mathcal{A}_{1p}\\
\vdots & \ddots & \vdots\\
\mathcal{A}_{p1} & \cdots & \mathcal{A}_{pp}
\end{bmatrix}.
$$

块耦合算子（block-coupled operator）定义其第 $i$ 个输出字段为

$$
y_i=\sum_{j=1}^{p}\mathcal{A}_{ij}x_j.
$$

每个子块 $\mathcal{A}_{ij}$ 可以是显式稠密矩阵、显式稀疏矩阵或矩阵自由算子。`BlockLayout` 只负责把 Krylov 所用的扁平向量映射到连续字段切片；`BlockOperator` 负责按上式累加子块作用。因此它保留字段边界，却不要求 Krylov 后端改用嵌套向量表示。

**Vlasov-Fokker-Planck 高阶矩—Maxwell 例子。** 可取

$$
x=\begin{bmatrix}m\\E\\B\end{bmatrix},
\qquad
\mathcal{A}=
\begin{bmatrix}
\mathcal{A}_{mm} & \mathcal{A}_{mE} & \mathcal{A}_{mB}\\
\mathcal{A}_{Em} & \mathcal{A}_{EE} & \mathcal{A}_{EB}\\
\mathcal{A}_{Bm} & \mathcal{A}_{BE} & \mathcal{A}_{BB}
\end{bmatrix},
$$

其中 $m$ 表示所选阶数的矩变量。对角块描述各字段自身的输运、碰撞或电磁离散；非对角块描述电流、洛伦兹力、旋度和本构耦合。实际块划分必须来自物理离散和边界条件，不能仅按变量名称猜测。

块结构本身不证明对称性、正定性或可逆性，也不自动构成预条件器。整体 Krylov 路线仍按整个算子的数学契约判定。其价值在于为后续字段分裂、块三角近似和 Schur 补预条件器提供可审计的字段边界：例如消去 $m$ 后，场变量的形式 Schur 补为

$$
\mathcal{S}_{EB}=
\mathcal{A}_{EB,EB}-\mathcal{A}_{EB,m}\mathcal{A}_{mm}^{-1}\mathcal{A}_{m,EB}.
$$

该表达只说明目标结构，不授权当前实现显式形成 $\mathcal{A}_{mm}^{-1}$ 或 $\mathcal{S}_{EB}$。

<a id="scale-and-conditioning"></a>
### 3.3 规模、资源与条件信息

矩阵规模会改变路线的经济性，但不能单独决定路线。`small`、`medium` 和 `large` 应由 `RouteBudget` 的内存、时间、右端项数、稀疏模式和硬件能力共同定义，而不应在核心中写死单一维度阈值。

| 资源类别 | 典型判据 | 条件信息策略 | 初始数学路线 |
|---|---|---|---|
| 小型 | 稠密分解可在预算内完成 | 可按需计算 $κ_p(A)$；也可从直接法获得诊断 | 主元 LU；病态或秩亏风险时 QR/SVD |
| 中型 | 稀疏直接分解可能可行，但 fill-in 或多右端决定代价 | 优先复用分解的倒数条件数与误差界；必要时进行低成本估计 | 多右端或高复用时直接法；否则预条件 Krylov |
| 大型 | 稀疏、矩阵自由、分布式或直接分解超预算 | 默认不做完整条件数计算；只接受先验、按需估计或迭代诊断 | 结构合格的 Krylov 加预条件器 |

条件数是可选输入，而不是必须由路由器重新计算。其数学语义必须完整保存：

```julia
ConditioningInfo(
    estimate = nothing,               # 条件数或倒数条件数的数值
    metric = :kappa_2,                # :kappa_1, :kappa_2, :kappa_inf, :rcond
    operator = :original,             # :original, :left_preconditioned, :right_preconditioned
    matrix_version = nothing,
    source = :caller_estimate,
    reliable = true,
    evidence = :external,             # :exact, :estimated, :external, :qualitative
)
```

例如，原始算子 $A$ 的条件数很大，而预条件算子 $A M^{-1}$ 的条件数可能适合 Krylov 迭代。路由性能判断应优先使用与候选路线一致的算子信息。当前实现要求 `matrix_version` 与 `AdaptiveLinearProblem.matrix_version` 一致，并要求调用方标记 `reliable=true`；过期、未知范数、未知预条件形式或不可信的信息会被记录为 `:stale`、`:invalid` 或 `:untrusted`，不会改变路线排序。已知条件数也不能绕过 Hermitian、正定、维度和伴随等算法资格。

`ConditioningPolicy()` 的默认 `estimation=:none`，因此只验证外部信息，不计算条件数。只有调用方显式设置 `estimation=:cheap` 或 `:full`，并同时给出正的 `DiagnosticBudget.max_seconds` 与 `max_matrix_dimension`，`diagnose` 才允许在小型显式稠密矩阵上估计条件信息。`cheap` 记录一范数倒数条件数，`full` 记录二范数条件数；矩阵自由、稀疏或超出维度预算的情况返回受限原因，而不是隐式稠密化。

### 3.2 条件数的计算与估计

条件数定义为

$$
\kappa_p(A)=\lVert A\rVert_p\lVert A^{-1}\rVert_p,\qquad
\operatorname{rcond}_p(A)=\kappa_p(A)^{-1}.
$$

对于小型显式矩阵，Julia 的 `cond(A, p)` 可计算 $p=1,2,\infty$ 的条件数；$p=2$ 路线依赖奇异值信息，适合诊断而不应成为大型问题的默认前置计算[[1](#ref-1)]。

对于已经选择直接分解的中型问题，条件诊断应尽量作为分解的副产物取得。LAPACK 的 `gesvx!` 可返回倒数条件数 `rcond`、前向误差界 `ferr`、后向误差界 `berr` 和主元增长信息；`gecon!` 可在已有 LU 分解上估计一范数或无穷范数的倒数条件数[[1](#ref-1)]。稀疏直接后端只有在实际暴露此类诊断时才记录它，不能假设每个后端都能低成本给出可靠估计。

对于大型 SPD 系统，可从 Lanczos 或 CG 过程中得到 Ritz 值并估计谱区间，从而给出条件数趋势；这属于运行时诊断，不是严格证书。对于一般矩阵，可通过一范数逆估计器配合 $A$ 与 $A^\ast$ 的求解估计 $\lVert A^{-1}\rVert_1$，但它本身需要额外求解，因而只应在诊断预算明确允许时使用。矩阵自由的一般非正规系统不存在同时廉价、可靠且通用的完整条件数前置估计。

当前实现以 `spectral_estimation=:lanczos` 将这一原则落实为受限路径：仅当 Hermitian 性已由 `Certified` 或 `Proved` 证据支持、矩阵为方阵、时间预算为正且允许正的算子调用数时才运行。算法只调用 $y=Ax$，因而可用于显式稀疏矩阵和矩阵自由算子；`SpectralInfo` 保存 Ritz 区间、步数、矩阵版本和来源。对已认证正定系统，极小的估计下端点只生成 `:near_singular` 诊断提示，不是秩亏证明，也不能自行触发或禁止数学路线。细节见[大规模条件与谱诊断](docs/SPECTRAL_DIAGNOSTICS.md)。

当条件信息未知时，`unknown` 不等于 `ill-conditioned`。路由器采用以下保守策略：

- 小型问题先选有主元的直接法，验收后向误差并记录可取得的条件诊断；若秩亏、误差或主元增长可疑，回退到 QR 或 SVD。
- 中型稀疏问题根据 fill-in 风险、右端项复用和内存预算选择直接法或预条件 Krylov；用短预算迭代的残差下降率决定是否重建预条件器或切换路线。
- 大型或矩阵自由问题按严格结构资格选择 Krylov 法与预条件器，不先计算完整条件数；持续监控残差、停滞、breakdown、预条件器代价和内存。

诊断输出必须分开保存条件、迭代和预条件器状态。`NumericalDiagnosis.conditioning.state` 可为 `:moderate`、`:ill_conditioned`、`:near_rank_deficient`、`:stale`、`:untrusted` 或 `:unknown`；`iteration_state` 可为 `:converged`、`:stagnated`、`:breakdown`、`:nonconverged` 或 `:not_observed`；`preconditioner_state` 只将停滞或 breakdown 标记为 `:suspected_failure`，不把相关性误写成因果结论。这样近秩亏、迭代停滞、预条件器可疑失效和实际未收敛具有独立可审计字段。

这与 `LinearSolve.jl` 接受调用方提供算子条件假设的思想一致，但本项目保留数值、范数、作用算子和证据等级，避免把单一的“良态/病态”标签误用于不匹配的路线[[6](#ref-6)]。

## 4. 受控自适应

自适应由三步组成：

1. 结构与资格筛选：排除数学上不合法的路线。
2. 主路线执行：在迭代、时间和内存预算内监控收敛。
3. 有限回退：仅在可解释事件后切换，例如直接分解失败、CG 检出非正曲率、GMRES 停滞、预条件器失效或预算耗尽。

一个问题默认不应尝试超过主路线和两条回退路线。切换路线时必须保存每次尝试的参数、证据、停止原因与残差历史。

<a id="fixed-and-variable-preconditioners"></a>
### 4.1 固定与可变预条件器

设右预条件系统为

$$
A M^{-1} y=b,\qquad x=M^{-1}y.
$$

GMRES 假定同一次线性求解的每一步均使用固定的线性预条件器 $M^{-1}$，并据此构造 Krylov 子空间。若第 $k$ 步实际使用 $M_k^{-1}$，固定预条件的 Arnoldi 关系不再成立；普通 GMRES 仍可能运行，但其标准残差最小化解释与收敛判据不再有保证。

FGMRES 允许每步使用不同的预条件器。它保留

$$
z_k=M_k^{-1}v_k,
$$

并在由各个 $z_k$ 构成的搜索空间中执行残差最小化，因此适用于可变或非线性预条件器[[5](#ref-5)]。

在一次 `solve(A, b)` 内出现下列情形时，规划器必须选择 FGMRES：

- ILU 参数、稀疏因子或块近似在 Krylov 迭代期间重建。
- AMG 层次、smoother 或粗网格求解器在 Krylov 迭代期间调整。
- 预条件器由内层迭代法实现，且内层迭代次数或容差随外层残差变化。
- 预条件器依赖当前外层 Krylov 向量，因而不是固定线性算子。

预条件器只在不同时间步、不同 Newton 步或不同的 `solve(A, b)` 调用之间更新时，每次单独求解内部仍可将其视为固定；此时普通 GMRES 仍具有资格。能力模型应显式记录 `fixed_within_solve` 与 `linear_within_solve` 两个布尔性质；任一性质未知或为假时，保守地选择 FGMRES。

## 5. 正确性与验收

数值成功至少要求报告：

- 相对残差 $\lVert b-Ax\rVert/\lVert b\rVert$。
- 后向误差 $\lVert b-Ax\rVert/(\lVert A\rVert\lVert x\rVert+\lVert b\rVert)$。
- 算法状态：收敛、停滞、breakdown、资格拒绝或预算耗尽。
- 迭代数、时间、内存和预条件器构造/应用代价。
- 路线证书及每项数学证据。

若调用方提供额外的数学约束，例如零空间投影、约束残差或加权内积，验收器将其作为额外条件，而非物理专用逻辑。

### 5.1 可观测性、历史与性能预算

求解器必须产出足以比较不同路线的结构化记录，才能在后续相似问题上提出建议。记录应保存矩阵特征、数学证据、路线、预条件器、停止原因和效果，而不是默认保存原始矩阵、右端项、解向量或完整 Krylov 基。

| 遥测级别 | 必要输出 | 额外开销 | 默认用途 |
|---|---|---|---|
| `off` | 解与瞬时 `SolveStatus`；不生成证书、不写入历史 | 无额外遥测开销 | 内层线性求解默认值 |
| `basic` | 调用方显式请求的固定字段，如路线、迭代数或最终残差；不构造指纹、不写入历史 | 通常为常数级；请求真实残差时至多增加一次算子应用 | 单次求解的定制输出 |
| `fingerprint` | 按标签与事件选择的字段，如矩阵指纹、预条件器代价、回退原因和环境摘要 | 通常很低；仅在触发事件时构造并输出 | 历史建议与相似问题统计 |
| `trace` | `fingerprint` 加采样的残差、时间和内存轨迹 | 与采样次数成比例；重新计算真实残差可能增加算子应用 | 调参与失败分析 |
| `diagnostic` | `trace` 加谱、条件数、主元增长或额外验证 | 可能显著，必须显式启用并受预算限制 | 数值研究与异常诊断 |

默认 `TelemetryPolicy` 为 `level=:off` 且 `emit_on=()`。它只返回解和求解状态，不写入 `HistoryStore`，不保存轨迹，也不为遥测重新计算残差。`basic`、`fingerprint`、`trace` 与 `diagnostic` 不应强制每一步重新计算残差。若后端已经维护递推残差，则直接记录该值；若用户要求真实残差，验收器才执行 $b-Ax$ 的额外计算。`trace` 默认只在预设间隔采样，不保存全部迭代向量。条件数估计、谱估计和高频内存采样属于 `diagnostic`，不进入默认路径。

`OutputRequest` 与 `FingerprintProfile` 是两个不同对象。前者指定本次调用无条件需要的字段，适合 `basic`；后者定义“标签组合 -> 输出字段”的映射，并与 `emit_on` 配合，适合 `fingerprint` 及以上级别。因而，“需要特定输出”不必然要求指纹模式；只有输出字段取决于问题标签，或记录将用于相似问题历史建议时，才使用指纹模式。

历史记录由标签指纹和事件共同控制。标签只能使用已经可得或由调用方提供的低成本信息；如果取得某个特征需要扫描矩阵、额外分解或额外算子应用，该标签保持 `unknown`，不得为收集历史而增加计算。

| 指纹标签 | 典型取值 | 事件发生时的对应输出 |
|---|---|---|
| 表示 | `dense`、`sparse`、`operator`、`block` | 维度、非零元数或算子接口能力 |
| 规模带 | `small`、`medium`、`large` | 预算类别与直接法可行性判断 |
| 数学结构 | `spd`、`hermitian_indefinite`、`general`、`rectangular` | 资格证据与被排除的方法 |
| 条件信息 | `unknown`、`hinted`、`estimated`、`very_ill` | 条件数对象、范数、版本与估计成本 |
| 预条件语义 | `none`、`fixed`、`variable` | 构造/应用代价、缓存命中与 GMRES/FGMRES 选择 |
| 执行形态 | `one_rhs`、`multi_rhs`、`matrix_free`、`distributed` | 复用信息、乘法次数或分布式后端摘要 |
| 问题族 | 调用方的 `family_key` | 仅用于相似历史检索，不替代数学资格 |

推荐的历史策略是 `level=:fingerprint`，并设置 `emit_on=(:new_fingerprint, :failure, :fallback, :periodic)`。这会为首次出现的标签组合、失败和回退保留完整摘要，而不是为每个内层线性系统写入一条记录。只有研究特定迭代行为时，才显式使用 `emit_on=(:every_solve,)` 或 `level=:trace`。

每条由事件触发、可用于历史学习的 `SolveRecord` 至少包含：

- 问题指纹：维度、非零元数、密度、标量类型、显式/矩阵自由表示、已认证结构、块布局摘要和调用方提供的 `family_key`。
- 条件信息：`ConditioningInfo`、条件估计的对象、范数、证据与矩阵版本。
- 路线信息：候选、最终方法、后端、预条件器、缓存命中、容差和预算。
- 结果信息：状态、迭代数、时间、内存摘要、最终残差、后向误差、回退原因。
- 环境信息：Julia 版本、包版本、设备类别和线程配置。

历史学习只能在资格筛选之后排序已合法路线。它可以推荐“对同一 `family_key` 和相似指纹，某预条件器曾有效”，但不能因历史成功而允许未认证 SPD 系统使用 CG，或允许缺少伴随的算子使用最小二乘算法。初始实现采用有界的内存记录库与确定性排序；持久化存储、统计模型和受控探索属于后续扩展。

当前实现将 `FingerprintProfile` 作为字段级开关：未启用的标签保持 `nothing`，不会参与相似记录检索。`HistoryStore` 仅在调用方传入实例且事件匹配时写入；`route_advice` 只对已通过资格门的自动候选路线提出重排建议。`save_history` 与 `load_history!` 必须由调用方显式调用，并且二进制历史只能从受信任本机路径读取。`trace` 与 `diagnostic` 只截取后端已有残差历史，受 `TelemetryBudget` 的样本上限约束；它们不会因遥测而额外计算真实残差、谱或条件数。完整接口见[可观测性与历史自适应](docs/OBSERVABILITY.md)。

资源预算在数学资格筛选之后、数值内核执行之前生效。时间和迭代上限可收紧 Krylov 控制；内存预算只在当前后端能够给出保守下界时通过，未知 fill-in 或不支持的设备形态必须返回可解释的预算拒绝。后端登记应区分“已知名称”“已安装依赖”和“已验证适配器”：只有已验证适配器才可执行。当前实现只有串行 CPU 的标准库直接法与 `Krylov.jl` 后端；GPU、MPI、PETSc、LinearSolve.jl 和 IterativeSolvers.jl 均不得因名称匹配而被假定可用。完整能力表与适配器准入条件见[计算资源与后端](docs/BACKENDS.md)。

## 6. 生态定位

`LinearSolve.jl` 提供统一线性问题接口、算法适配与缓存；本项目将它作为可选执行桥，而不是将其默认启发式视为数学证明[[3](#ref-3)]。`AlgebraicMultigrid.jl` 适合以扩展形式接入已认证的稀疏 SPD 路线[[4](#ref-4)]。PETSc、GPU 和分布式内存后端均属于后续扩展，不进入数学核心。

## 参考文献列表

[1]. <a id="ref-1"></a> Julia Language. 2026. [*Linear Algebra*](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/). Julia Documentation.<br>

[2]. <a id="ref-2"></a> JuliaSmoothOptimizers. 2026. [*Matrix-free operators*](https://github.com/JuliaSmoothOptimizers/Krylov.jl/blob/main/docs/src/matrix_free.md). Krylov.jl Documentation.<br>

[3]. <a id="ref-3"></a> SciML. 2026. [*Algorithm Selection Guide*](https://docs.sciml.ai/LinearSolve/dev/basics/algorithm_selection/). LinearSolve.jl Documentation.<br>

[4]. <a id="ref-4"></a> JuliaLinearAlgebra. 2026. [*Algebraic Multigrid in Julia*](https://github.com/JuliaLinearAlgebra/AlgebraicMultigrid.jl). GitHub repository.<br>

[5]. <a id="ref-5"></a> JuliaSmoothOptimizers. 2026. [*Reference*](https://jso.dev/Krylov.jl/dev/interfaces/reference/). Krylov.jl Documentation.<br>

[6]. <a id="ref-6"></a> SciML. 2026. [*Linear Solve Operator Assumptions*](https://docs.sciml.ai/LinearSolve/v2.30/basics/OperatorAssumptions/). LinearSolve.jl Documentation.<br>
