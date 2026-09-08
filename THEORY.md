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

对病态系统，路由器不应常规计算完整条件数。完整条件数通常昂贵，且对非正规矩阵的 Krylov 收敛预测能力有限。更可靠的运行时指标是后向误差、预条件后残差下降率、停滞、breakdown、内存占用和用户规定的预算。

## 4. 受控自适应

自适应由三步组成：

1. 结构与资格筛选：排除数学上不合法的路线。
2. 主路线执行：在迭代、时间和内存预算内监控收敛。
3. 有限回退：仅在可解释事件后切换，例如直接分解失败、CG 检出非正曲率、GMRES 停滞、预条件器失效或预算耗尽。

一个问题默认不应尝试超过主路线和两条回退路线。切换路线时必须保存每次尝试的参数、证据、停止原因与残差历史。

## 5. 正确性与验收

数值成功至少要求报告：

- 相对残差 $\lVert b-Ax\rVert/\lVert b\rVert$。
- 后向误差 $\lVert b-Ax\rVert/(\lVert A\rVert\lVert x\rVert+\lVert b\rVert)$。
- 算法状态：收敛、停滞、breakdown、资格拒绝或预算耗尽。
- 迭代数、时间、内存和预条件器构造/应用代价。
- 路线证书及每项数学证据。

若调用方提供额外的数学约束，例如零空间投影、约束残差或加权内积，验收器将其作为额外条件，而非物理专用逻辑。

## 6. 生态定位

`LinearSolve.jl` 提供统一线性问题接口、算法适配与缓存；本项目将它作为可选执行桥，而不是将其默认启发式视为数学证明[[3](#ref-3)]。`AlgebraicMultigrid.jl` 适合以扩展形式接入已认证的稀疏 SPD 路线[[4](#ref-4)]。PETSc、GPU 和分布式内存后端均属于后续扩展，不进入数学核心。

## 参考文献列表

[1]. <a id="ref-1"></a> Julia Language. 2026. [*Linear Algebra*](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/). Julia Documentation.<br>

[2]. <a id="ref-2"></a> JuliaSmoothOptimizers. 2026. [*Matrix-free operators*](https://github.com/JuliaSmoothOptimizers/Krylov.jl/blob/main/docs/src/matrix_free.md). Krylov.jl Documentation.<br>

[3]. <a id="ref-3"></a> SciML. 2026. [*Algorithm Selection Guide*](https://docs.sciml.ai/LinearSolve/dev/basics/algorithm_selection/). LinearSolve.jl Documentation.<br>

[4]. <a id="ref-4"></a> JuliaLinearAlgebra. 2026. [*Algebraic Multigrid in Julia*](https://github.com/JuliaLinearAlgebra/AlgebraicMultigrid.jl). GitHub repository.<br>
