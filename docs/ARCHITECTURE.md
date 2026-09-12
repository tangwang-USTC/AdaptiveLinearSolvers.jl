# AdaptiveOptimization.jl 自适应优化求解器架构设计

## 版本与状态

| 版本 | 日期 | 作者 | 变更 |
|---|---|---|---|
| v0.1 | 2026-09-12 | 王彦鹏 | 初始架构设计 |

## 术语与缩略词

| 缩略词 | 首次定义 |
|---|---|
| LP | Linear Programming，线性规划 |
| QP | Quadratic Programming，二次规划 |
| NLP | Nonlinear Programming，非线性规划 |
| SDP | Semidefinite Programming，半定规划 |
| SOCP | Second-Order Cone Programming，二阶锥规划 |
| MILP | Mixed-Integer Linear Programming，混合整数线性规划 |
| KKT | Karush-Kuhn-Tucker conditions，卡罗需-库恩-塔克条件 |
| SQP | Sequential Quadratic Programming，序列二次规划 |
| IPM | Interior Point Method，内点法 |
| CG | Conjugate Gradient，共轭梯度法 |
| LBFGS | Limited-memory BFGS，有限内存 BFGS |
| TR | Trust Region，信赖域 |
| AD | Automatic Differentiation，自动微分 |
| MOI | MathOptInterface，数学优化接口 |
| SPD | Symmetric Positive Definite，对称正定 |

---

## 1. 愿景与范围

### 1.1 项目定位

构建一个统一的 Julia 优化求解器框架 `AdaptiveOptimization.jl`。它以统一接口接收优化问题，根据问题结构、数学资格条件、资源预算和可选诊断信息，在**无约束方法、约束优化方法、线性/二次规划求解器、非线性规划求解器、最小二乘求解器**之间进行可解释的自适应路由；用户也可在任意层级指定路线，其余层级仍采用默认或自适应策略。

### 1.2 与 AdaptiveLinearSolvers.jl 的关系

```
AdaptiveOptimization.jl          ← 本架构描述的优化求解器框架
  ├── core/                       ← 问题表示、路由、证书
  ├── algorithms/                 ← 优化算法家族
  │     └── ... (梯度下降、Newton、SQP、IPM ...)
  └── linearsolve/                ← Auto 模式委托给 AdaptiveLinearSolvers.jl
        └── AdaptiveLinearSolvers.jl   ← 作为内层线性系统求解器
```

优化算法内部需要反复求解线性系统（KKT 系统、Newton 方程、最小二乘正规方程），这部分全部委托给 `AdaptiveLinearSolvers.jl` 或其等价内核。因此 `AdaptiveLinearSolvers.jl` 的**全部数学资格门、路由策略、预条件器契约、资源预算和证书体系**在本框架中直接继承。

### 1.3 服务目标

- **Vlasov-Fokker-Planck 方程** 隐式时间步中的无约束/约束优化子问题
- **高阶矩方程组** 的约束投影与最小化
- **Maxwell 方程** 逆问题与参数估计
- **King mixture 逼近** 的似然最大化与正则化
- **等离子体平衡重建** 的约束优化

---

## 2. 设计原则

继承 `AdaptiveLinearSolvers.jl` 的全部五条设计原则并扩展：

1. **数学核心不依赖物理领域语义。** 优化问题和算法定义只使用数学对象：目标函数、约束函数、变量边界；不嵌入物理单位、场方程或离散格式。

2. **路由资格先于性能启发式。** 算法选择必须有数学资格证据（凸性、光滑性、约束品格），而不是仅凭问题规模或历史性能。无资格证据时默认走保守路线。

3. **用户可以锁定任意决策层，其余层保持自动。** 继续使用 `Auto / Prefer / Lock / Forbid` 四值策略，可覆盖家族、具体算法、容差、预条件器、后端等各层。

4. **后端是可选能力，不能污染核心 API。** Jul ia 原生算法、外部求解器、GPU/分布式后端均通过能力模型注册，框架核心不假设任何特定后端的可用性。

5. **每次结果必须带可机器读取和可人工审查的路线证书。** 证书记录候选算法、资格证据、尝试顺序、停止原因、最终残差/最优性度量，以及诊断信息。

**新增原则：**

6. **求解器输出必须可验证。** 框架提供统一的最优性验证器（KKT 残差、对偶间隙、约束违反度），不依赖后端自身的成功标志。

7. **导数信息是契约的一部分。** 自动微分、解析梯度、有限差分和矩阵自由作用均通过 `DerivativeContract` 注册，算法根据可用导数信息选择自身版本。

---

## 3. 问题表示层

### 3.1 统一优化问题

```julia
struct OptimizationProblem{TV, TO, TC, TH, TM}
    # 变量
    n_vars::Int
    has_integer::Bool
    integer_mask::Union{Nothing, BitVector}
    lower_bounds::Union{Nothing, TV}
    upper_bounds::Union{Nothing, TV}

    # 目标
    objective::ObjectiveFunction
    objective_type::Symbol  # :linear, :quadratic, :nonlinear, :least_squares

    # 约束
    constraints::ConstraintSet
    n_eq::Int
    n_ineq::Int

    # 数学契约
    contract::OptimizationContract

    # 导数契约
    derivatives::DerivativeContract

    # 提示与元数据
    hints::OptimizationHints
    label::Symbol
    matrix_version::Any
end
```

### 3.2 目标函数

```julia
@enum ObjectiveType begin
    Linear
    Quadratic
    NonlinearLeastSquares
    NonlinearGeneral
end

struct ObjectiveFunction
    f::Function                    # f(x) → scalar
    ∇f!::Union{Nothing, Function}  # ∇f!(g, x) → nothing, g overwritten
    Hv!::Union{Nothing, Function}  # Hv!(Hv, x, v) → nothing (Hessian-vector product)

    # 结构化目标
    linear_terms::Union{Nothing, AbstractVector}    # c  for cᵀx
    quadratic_operator::Union{Nothing, AbstractMatrix}  # Q for ½xᵀQx
    residual!::Union{Nothing, Function}  # r!(r, x) for ½‖r(x)‖²
    jacobian!::Union{Nothing, Function}  # J!(J, x) for residual Jacobian
end
```

### 3.3 约束系统

```julia
struct ConstraintSet
    # 线性约束: A x ◇ b
    linear_eq::Union{Nothing, AbstractMatrix}
    linear_eq_rhs::Union{Nothing, AbstractVector}
    linear_ineq::Union{Nothing, AbstractMatrix}
    linear_ineq_lhs::Union{Nothing, AbstractVector}
    linear_ineq_rhs::Union{Nothing, AbstractVector}

    # 非线性约束: c(x) ◇ 0
    nonlinear_eq::Union{Nothing, Function}     # c_eq(x) → vector
    ∇nonlinear_eq::Union{Nothing, Function}    # ∇c_eq!(J, x) → nothing
    nonlinear_ineq::Union{Nothing, Function}   # c_ineq(x) → vector
    ∇nonlinear_ineq::Union{Nothing, Function}  # ∇c_ineq!(J, x) → nothing

    # 锥约束
    cone::Symbol  # :none, :second_order, :semidefinite, :exponential
    cone_dimensions::Vector{Int}
end
```

### 3.4 数学契约（扩展）

在 `AdaptiveLinearSolvers.MathematicalContract` 的基础上扩展优化领域的性质：

```julia
Base.@kwdef struct OptimizationContract
    # 目标性质
    objective_convex::PropertyEvidence = PropertyEvidence()
    objective_strongly_convex::PropertyEvidence = PropertyEvidence()
    objective_quadratic::PropertyEvidence = PropertyEvidence()
    objective_self_concordant::PropertyEvidence = PropertyEvidence()
    objective_lipschitz_gradient::PropertyEvidence = PropertyEvidence()
    lipschitz_constant::Union{Nothing, Float64} = nothing

    # 约束性质
    constraints_linear::PropertyEvidence = PropertyEvidence()
    constraints_convex::PropertyEvidence = PropertyEvidence()
    constraints_slater::PropertyEvidence = PropertyEvidence()  # Slater's condition
    constraint_qualification::Symbol = :unknown  # :LICQ, :MFCQ, :Slater, :unknown

    # 全局性质
    problem_class::Symbol = :unknown  # :LP, :QP, :NLP, :least_squares, etc.
    partially_separable::PropertyEvidence = PropertyEvidence()
end
```

### 3.5 导数契约

```julia
@enum ADMode begin
    ForwardAD
    ReverseAD
    FiniteDifference
    Analytical
    MatrixFree
end

Base.@kwdef struct DerivativeContract
    gradient_mode::Union{Nothing, ADMode} = nothing
    hessian_mode::Union{Nothing, ADMode} = nothing
    jacobian_mode::Union{Nothing, ADMode} = nothing
    hv_product_available::Bool = false          # Hessian-vector product
    jv_product_available::Bool = false          # Jacobian-vector product
    jtv_product_available::Bool = false         # Jacobianᵀ-vector product

    # 精度保证
    gradient_precise::Bool = false
    hessian_exact::Bool = false
end
```

### 3.6 问题提示

```julia
Base.@kwdef struct OptimizationHints
    # 规模提示
    n_vars_hint::Symbol = :unknown    # :small (<100), :medium (100-10K), :large (>10K)
    n_constraints_hint::Symbol = :unknown
    density_hint::Symbol = :unknown   # :dense, :sparse, :matrix_free

    # 数值提示
    expected_conditioning::Symbol = :unknown
    scaling_issue::Bool = false
    degenerate::PropertyEvidence = PropertyEvidence()

    # 性能提示
    multiple_rhs::Bool = false
    hot_start_available::Bool = false
    warm_start_point::Union{Nothing, AbstractVector} = nothing
end
```

---

## 4. 问题分析与分类层

### 4.1 自动结构检测

问题创建后，框架运行一组低成本的结构检测器，用结果补充调用方提供的 `OptimizationContract`：

```julia
struct ProblemAnalysis
    # 确定性分类（低计算成本）
    # 返回 ClassificationResult{conclusion, evidence_level}
    problem_class::ClassificationResult  # :LP, :QP, :NLP, :least_squares...
    has_bounds::Bool
    has_eq_constraints::Bool
    has_ineq_constraints::Bool
    has_integer_vars::Bool

    # 推断性质（低到中等成本，以 contract 形式存储）
    inferred_contract::OptimizationContract

    # 规模指标
    effective_n::Int         # 实际决策变量数
    total_constraints::Int   # 约束总数（含边界）
    nnz_objective::Union{Nothing, Int}
    nnz_jacobian::Union{Nothing, Int}

    # 诊断
    notes::Vector{String}
end
```

**检测优先级**（成本从低到高）：

| 检测项 | 成本 | 依赖 | 产出 |
|---|---|---|---|
| 变量结构 | O(n) | 边界、整数标记 | has_bounds, has_integer |
| 线性约束识别 | O(nnz) | A 矩阵 | 约束线性/非线性分类 |
| 目标分类 | O(1) 或 O(nnz) | 目标函数结构 | LP/QP/NLP/最小二乘 |
| 凸性推断 | 1 次 Hessian 乘法 | Hv! 或显式 Hessian | 强凸/凸/非凸 |
| 稀疏模式 | O(nnz) | 显式矩阵 | 密度提示 |
| 条件估计 | 按预算 | 显式矩阵 | 数值提示 |

### 4.2 问题分类体系

```
OptimizationProblem
  ├── Continuous
  │     ├── Unconstrained
  │     │     ├── Smooth (C²)
  │     │     │     ├── Convex (包括强凸)
  │     │     │     ├── Nonconvex (一般)
  │     │     │     └── Least Squares
  │     │     └── Nonsmooth
  │     │           ├── L1-regularized
  │     │           └── Composite
  │     └── Constrained
  │           ├── Linearly Constrained
  │           │     ├── LP
  │           │     ├── QP (凸/非凸)
  │           │     └── Bounded (仅边界约束)
  │           ├── Nonlinearly Constrained
  │           │     ├── 凸 NLP
  │           │     └── 一般 NLP
  │           └── Conic
  │                 ├── SOCP
  │                 └── SDP
  └── Mixed-Integer
        ├── MILP
        ├── MIQP
        └── MINLP
```

---

## 5. 求解器注册与能力模型

### 5.1 求解器能力描述

```julia
Base.@kwdef struct SolverCapability
    # 标识
    name::Symbol
    package::Symbol
    version::VersionNumber

    # 解决的问题类
    problem_classes::Vector{Symbol}   # :LP, :QP, :NLP, :least_squares...

    # 变量类型
    supports_integer::Bool
    supports_bounded::Bool

    # 约束类型
    supports_eq_constraints::Bool
    supports_ineq_constraints::Bool
    supports_conic_constraints::Bool

    # 导数需求
    required_derivative::Symbol    # :none, :gradient, :hessian, :jacobian
    supports_matrix_free::Bool

    # 执行模式
    execution_modes::Vector{Symbol}  # :serial_cpu, :gpu, :distributed
    large_scale::Bool

    # 实现状态
    implemented::Bool     # true = 当前可调用
    available::Bool       # true = 已安装依赖

    # 数值特性
    exact_solver::Bool    # true = 返回精确解（LP/QP 等）
    global_guarantee::Symbol  # :global, :local, :none
end
```

### 5.2 算法注册表（初始）

| 求解器名称 | 问题类 | 整数 | 约束 | 导数需求 | 实现状态 |
|---|---|---|---|---|---|
| `steepest_descent` | 无约束光滑 | 否 | 无 | 梯度 | 待实现 |
| `lbfgs` | 无约束光滑 | 否 | 无 | 梯度 | 待实现 |
| `newton` | 无约束光滑 C² | 否 | 无 | Hessian | 待实现 |
| `trust_region_newton` | 无约束光滑 C² | 否 | 无 | Hessian | 待实现 |
| `cg_descent` | 无约束凸光滑 | 否 | 无 | 梯度 | 待实现 |
| `nelder_mead` | 无约束非光滑 | 否 | 无 | 无 | 待实现 |
| `lbfgs_b` | 有界 | 否 | 边界 | 梯度 | 待实现 |
| `ipm_lp` | LP | 否 | 线性 | 无（已结构化） | 待实现 |
| `ipm_qp` | QP | 否 | 线性 | 无（已结构化） | 待实现 |
| `sqp` | NLP | 否 | 非线性 | 梯度+Jacobian | 待实现 |
| `ipopt` | NLP | 否 | 非线性 | 梯度+Jacobian+Hessian | 外部 |
| `osqp` | QP | 否 | 线性 | 无（已结构化） | 外部 |
| `gurobi` | LP/QP/MILP/MIQP | 是 | 线性 | 无（已结构化） | 外部 |
| `hi gh s` | LP/MILP | 是 | 线性 | 无（已结构化） | 外部 |
| `least_squares_lm` | 最小二乘 | 否 | 无 | Jacobian | 待实现 |
| `least_squares_trl` | 最小二乘 | 否 | 无 | Jacobian | 待实现 |

### 5.3 求解器注册

```julia
# 注册一个求解器到全局注册表
register_solver!(capability::SolverCapability, solve_fn::Function)

# 查询可求解某问题的求解器
function matching_solvers(problem::OptimizationProblem,
                          analysis::ProblemAnalysis) -> Vector{SolverCapability}
    # 1. 按问题类过滤
    # 2. 按变量/约束类型过滤
    # 3. 按导数可用性过滤
    # 4. 按可用状态过滤
    # 5. 按实现状态过滤
end
```

---

## 6. 自适应路由引擎

### 6.1 策略体系（继承与扩展）

```julia
Base.@kwdef struct OptimizationRoutePolicy
    # 求解器家族
    family::RouteChoice = Auto()       # :unconstrained, :constrained, :least_squares

    # 具体算法（家族内选择）
    algorithm::RouteChoice = Auto()    # :lbfgs, :newton, :sqp, :ipm_lp...

    # 预条件器（用于内层线性系统）
    preconditioner::RouteChoice = Auto()

    # 线性求解器（委托给 AdaptiveLinearSolvers）
    linear_solver::RouteChoice = Auto()

    # 后端
    backend::RouteChoice = Auto()      # :native, :ipopt, :gurobi, :osqp...

    # 自动微分
    ad_mode::RouteChoice = Auto()      # :forward, :reverse, :finite_diff

    # 容差和预算
    tolerances::RouteChoice = Auto()

    # 回退
    fallback::RouteChoice = Auto()
end
```

### 6.2 路线规划

```julia
struct OptimizationRoutePlan
    candidate_algorithms::Vector{Symbol}
    eligibility::Vector{EligibilityDecision}
    planned_algorithms::Vector{Symbol}
    execution_order::Vector{Symbol}
    unavailable::Vector{Symbol}
    layer_decisions::Vector{LayerDecision}
    inner_linear_plan::Union{Nothing, RoutePlan}  # AdaptiveLinearSolvers 的路线规划
    diagnosis::Union{Nothing, OptimizationDiagnosis}
    history_advice::Union{Nothing, RouteAdvice}
end
```

**规划流程（与线性求解器一致的分阶段设计）：**

```
plan(problem, policy)
  ├── Analyze problem structure (分析层)
  ├── Match solvers (注册表 → 候选算法)
  ├── Apply policy layer decisions (Auto/Prefer/Lock/Forbid)
  ├── Check mathematical eligibility (资格门)
  ├── Rank eligible algorithms (排序)
  ├── Generate inner linear plans (委派给线性求解器的 plan)
  └── Return OptimizationRoutePlan
```

### 6.3 资格判定矩阵

| 算法 | 凸目标 | 光滑 C² | 约束类型 | 导数需求 | 特殊情况 |
|---|---|---|---|---|---|
| 最速下降 | 不需要 | 需要 | 无 | 梯度 | |
| LBFGS | 不需要 | 需要 | 无/边界 | 梯度 | 边界时需投影 |
| Newton | 不需要 | 需要 | 无 | Hessian | 非凸时需正则化 |
| 信赖域 Newton | 不需要 | 需要 | 无 | Hessian | 自然处理非凸 |
| CG_DESCENT | 需要 | 需要 | 无 | 梯度 | |
| IPM LP | LP 类 | N/A | 线性 | N/A | 隐含 KKT |
| IPM QP | 凸 QP | N/A | 线性 | N/A | 隐含 KKT |
| SQP | 不需要 | 需要 | 非线性 | 梯度+Jacobian | 每步解 QP 子问题 |
| Levenberg-Marquardt | 最小二乘 | 需要 | 无 | Jacobian | |

### 6.4 执行流程

```julia
function solve(problem::OptimizationProblem; kwargs...)
    # 1. 分析
    analysis = analyze(problem)

    # 2. 规划
    plan = optimize_plan(problem, analysis, policy)

    # 3. 执行（带回退）
    for algorithm in plan.execution_order
        result = try
            execute_algorithm(problem, algorithm, plan)
        catch error
            # 记录失败，继续下一候选
            continue
        end
        # 4. 验证
        if verify_optimality(problem, result)
            return result
        end
    end

    # 5. 全部失败 → 返回诊断
    return OptimizationSolution(; status=AllFailed, diagnosis=...)
end
```

### 6.5 内层线性求解器的委托

优化算法每步需要求解线性系统（Newton 方程、KKT 系统），全权委托给 `AdaptiveLinearSolvers.jl`：

```julia
function solve_newton_equation(H, g; kwargs...)
    # 将 Hessian H 和右端项 g 打包为 AdaptiveLinearProblem
    # 委托给 AdaptiveLinearSolvers.solve()
    # 返回解 + 路线证书（嵌入到优化路线的证书中）
end
```

| 优化算法 | 内层线性系统 | 典型结构 | 委托路线偏好 |
|---|---|---|---|
| Newton | $H ∇f(x) = -∇f(x)$ | SPD（凸）或对称不定 | 自动（Cholesky/LU → Krylov） |
| SQP | KKT 系统 | 鞍点矩阵（对称不定） | LU → MINRES → GMRES |
| IPM LP/QP | 增广 KKT 系统 | 对称不定，趋于病态 | 预条件 MINRES → 直接法 |
| LM | $(JᵀJ+μI)Δx = -Jᵀr$ | SPD | Cholesky → CG |
| 信赖域 | 子问题：$||Δx||≤Δ$ | 结构依赖 | 专有方法 |

---

## 7. 算法家族

### 7.1 无约束优化

```
Unconstrained
  ├── Line Search
  │     ├── Steepest Descent       —— 梯度 + Armijo/Wolfe 线搜索
  │     ├── CG (Fletcher-Reeves/PR+) —— 仅凸
  │     ├── LBFGS                  —— 有限内存，适合大规模
  │     ├── Newton                 —— 精确 Hessian / 拟 Newton
  │     └── Newton-CG              —— Hessian-free，CG 解 Newton 方程
  └── Trust Region
        ├── Cauchy Point           —— 简单折线法
        ├── Steihaug CG            —— 适用于凸问题
        └── Newton-TR              —— 精确/近似 Hessian + 信赖域子问题
```

**资格链示例（LBFGS）：**

1. 问题无约束或仅有边界约束 → `Eligible`
2. 目标光滑（梯度可用）→ `Eligible`
3. 约束非线性 → `Ineligible`（走 SQP 或 IPM）
4. 未提供梯度、AD 不可用且不允许有限差分 → `Ineligible`

### 7.2 约束优化

```
Constrained
  ├── Active Set
  │     ├── Simplex (LP)           —— Dantzig 单纯形法
  │     └── QP Active Set          —— 小规模凸 QP
  ├── Interior Point
  │     ├── IPM-LP                 —— 线性规划原对偶内点法
  │     ├── IPM-QP                 —— 凸二次规划内点法
  │     └── IPM-NLP                —— 非线性原对偶内点法
  ├── Sequential Quadratic Programming
  │     ├── SQP (Line Search)      —— 线搜索 SQP
  │     └── SQP (Trust Region)     —— 信赖域 SQP (不精确 SQP)
  ├── Augmented Lagrangian
  │     ├── LANCELOT               —— 增广 Lagrangian + 边界约束
  │     └── ADMM                   —— 交替方向乘子法
  └── Penalty Methods
        ├── Quadratic Penalty
        └── Exact Penalty (L1)
```

**资格链示例（SQP）：**

1. 一般约束 NLP → `Eligible`
2. 凸无约束 → `Ineligible`（给无约束方法）
3. LP/QP → `Ineligible`（给专用方法）
4. 导数：至少需要约束 Jacobian → `Available`
5. 内层 QP求解器可用 → `Eligible`
6. 非光滑目标 → `Ineligible`（给 bundle 或子梯度方法）

### 7.3 最小二乘

```
Least Squares
  ├── Linear LS
  │     ├── QR (稠密)
  │     ├── SVD (秩亏)
  │     ├── Normal Equations (Cholesky/CG) —— 仅条件良好时
  │     └── LSQR/LSMR (大型稀疏)
  └── Nonlinear LS
        ├── Gauss-Newton             —— 残差小时接近 Newton
        ├── Levenberg-Marquardt      —— 信赖域变体，处理秩亏
        └── NLF (NLLS)               —— 大残差时的 Newton 法
```

最小二乘的资格与路由通过 `AdaptiveLinearSolvers` 的线性最小二乘路线和本框架的非线性最小二乘路线共同覆盖。

### 7.4 混合整数

```
Mixed-Integer
  ├── MILP
  │     ├── Branch & Bound
  │     └── Branch & Cut
  ├── MIQP
  │     └── B&B / B&C
  └── MINLP
        ├── B&B (NLP relaxation)
        └── Outer Approximation
```

---

## 8. 内层线性求解器集成

### 8.1 接口契约

优化算法中产生的所有线性系统使用统一的内层求解接口：

```julia
# 优化框架内：委托给 AdaptiveLinearSolvers
function inner_solve(problem::OptimizationProblem, A, b;
                     structure::Symbol = :general,  # :spd, :symmetric, :general, :kkt
                     preconditioner = nothing,
                     linear_solver_policy = :auto)
    # 1. 构建 AdaptiveLinearProblem（自动推断数学契约）
    # 2. 根据 structure 设置 contract（Hermitian / SPD 等）
    # 3. 必要时附加 preconditioner
    # 4. 调用 AdaptiveLinearSolvers.solve()
    # 5. 返回解 + 路线证书
end
```

### 8.2 KKT 系统的特殊处理

KKT 系统是优化中最关键的内层线性系统：

$$
\begin{bmatrix}
H & A^T \\
A & 0
\end{bmatrix}
\begin{bmatrix}
\Delta x \\
\nu
\end{bmatrix}
=
\begin{bmatrix}
-g \\
-h
\end{bmatrix}
$$

对 KKT 系统，框架提供：

1. **正则化**：病态 KKT 系统的 Hessian 正则化 $H \leftarrow H + \delta I$
2. **预条件器**：`KKTPreconditioner` 族（块对角、约束残量、Schur 补近似）
3. **简化**：对纯等式约束，消去 $\nu$ 得到简化系统
4. **惰性分解**：同一 SQP/IPM 迭代中复用 KKT 分解

```julia
Base.@kwdef struct KKTSystem
    H::AbstractMatrix
    A::AbstractMatrix
    regularization::Float64 = 0.0
    reduced::Bool = false  # true = 已消去对偶变量
end

function build_kkt_preconditioner(kkt::KKTSystem, strategy::Symbol)
    # :block_diagonal  —— diag(H, -A Aᵀ) 近似逆
    # :constraint_residual —— 约束残量缩放
    # :schur_complement   —— 显式或近似 Schur 补
    # :user              —— 调用方提供
end
```

### 8.3 预条件器共享

优化迭代中，对外层算法而言，内层线性系统的预条件器并非独立选择，而是和外层求解状态密切相关：

| 外层算法 | 内层系统 | 预条件器策略 |
|---|---|---|
| Newton | $H_k d = -∇f_k$ | 不可用（仅当 H 已知正定） |
| Newton-CG | $H_k d = -∇f_k$ | Jacobi / 不完全 Cholesky |
| SQP | KKT 系统 | Schur 补 / 块对角 |
| IPM | 增广 KKT | 正规方程简化 / 约束预条件器 |
| LM | $(J_k^T J_k + μI) d = -J_k^T r_k$ | 不需要（系统本身正定） |

---

## 9. 验证与证书

### 9.1 最优性验证

所有路线返回的候选解必须经过统一验证：

```julia
struct OptimalityVerification
    # 原始可行性
    max_constraint_violation::Union{Nothing, Float64}
    bound_violation::Union{Nothing, Float64}

    # 对偶可行性
    dual_feasibility::Union{Nothing, Float64}    # ‖∇f - Aᵀy - z‖

    # 互补松弛
    complementarity::Union{Nothing, Float64}     # |zᵢ xᵢ|

    # KKT 综合
    kkt_residual::Union{Nothing, Float64}
    kkt_satisfied::Bool

    # 目标信息
    f_x::Union{Nothing, Float64}
    grad_norm::Union{Nothing, Float64}

    # 验证状态
    status::Symbol  # :optimal, :suboptimal, :infeasible, :unbounded, :not_verified
    notes::Vector{String}
end
```

### 9.2 优化路线证书

```julia
struct OptimizationCertificate
    # 问题摘要
    problem_class::Symbol
    n_vars::Int
    n_constraints::Int

    # 路线信息
    candidate_algorithms::Vector{Symbol}
    attempted::Vector{Symbol}
    selected::Union{Nothing, Symbol}
    fallback_reason::Union{Nothing, Symbol}

    # 解验证
    verification::OptimalityVerification

    # 内层线性求解证书（嵌套）
    inner_linear_certificates::Vector{RouteCertificate}

    # 迭代统计
    iterations::Int
    elapsed_seconds::Float64
    function_evaluations::Int
    gradient_evaluations::Int
    hessian_evaluations::Int
    linear_solves::Int

    # 诊断
    diagnosis::Union{Nothing, OptimizationDiagnosis}
    notes::Vector{String}
end
```

### 9.3 测试矩阵

| 问题 | 类型 | n | 约束 | 特征 | 期望算法 |
|---|---|---|---|---|---|
| Rosenbrock | 无约束 | 2 | 0 | 非凸峡谷 | LBFGS / Newton |
| Himmelblau | 无约束 | 2 | 0 | 多局部极小 | 信赖域 |
| LP 示例 | LP | 可变 | 线性 | 多面体 | IPM / Simplex |
| QP 示例 | 凸 QP | 可变 | 线性 | 正定 Hessian | IPM / OSQP |
| SQP 示例 | NLP | 可变 | 非线性 | 约束优化 | SQP |
| NLS 示例 | 最小二乘 | 可变 | 0 | 非线性残差 | LM / Gauss-Newton |

---

## 10. 可观测性与自适应

### 10.1 遥测体系

继承 `AdaptiveLinearSolvers` 的五级遥测层次：

| 级别 | 优化特定输出 | 开销 | 用途 |
|---|---|---|---|
| `off` | 解 + 状态 | 无 | 生产内层 |
| `basic` | 算法、迭代数、最终 f(x) | 常数 | 单次求解 |
| `fingerprint` | 问题指纹 + 算法选择 + 收敛历史 | 低 | 自适应学习 |
| `trace` | 函数值序列、梯度范数、步长 | 与迭代数成比例 | 收敛分析 |
| `diagnostic` | 谱、条件、KKT 诊断 | 高，需显式预算 | 研究 |

### 10.2 问题指纹（优化扩展）

```julia
Base.@kwdef struct OptimizationFingerprint
    # 基础（继承自 MatrixFingerprint）
    representation::Union{Nothing, Symbol}  # :explicit, :matrix_free
    size_band::Union{Nothing, Symbol}       # :small, :medium, :large

    # 优化特征
    problem_class::Union{Nothing, Symbol}   # :LP, :QP, :NLP, :least_squares
    convexity::Union{Nothing, Symbol}       # :strongly_convex, :convex, :nonconvex
    constraint_type::Union{Nothing, Symbol} # :unconstrained, :bound, :linear, :nonlinear
    has_integer::Union{Nothing, Bool}

    # 执行特征
    execution::Union{Nothing, Symbol}       # :serial_cpu, :gpu...
    ad_mode::Union{Nothing, Symbol}         # :forward, :reverse, :finite_diff
end
```

### 10.3 自适应学习

```julia
# 对相似指纹的优化问题，历史自适应提供：
struct OptimizationRouteAdvice
    matching_records::Int
    algorithm_scores::Dict{Symbol, Float64}
    recommended_algorithm::Union{Nothing, Symbol}
    recommended_linear_solver_policy::Union{Nothing, Symbol}
    recommended_preconditioner::Union{Nothing, Symbol}
    suggested_parameters::Dict{Symbol, Any}
    selection_reason::Symbol
end
```

---

## 11. 外部后端接口

### 11.1 统一封装层

```julia
# 使用 Julia 的 PackageCompiler 或直接调用模式：
abstract type ExternalSolverBackend end

struct IpoptBackend <: ExternalSolverBackend
    # Ipopt.jl 参数映射
    options::Dict{String, Any}
end

struct OSQPBackend <: ExternalSolverBackend
    # OSQP.jl 参数映射
    options::Dict{Symbol, Any}
end

# 统一求解接口
function solve_backend(backend::ExternalSolverBackend,
                       problem::OptimizationProblem,
                       options::Dict) -> OptimizationSolution
end
```

### 11.2 桥接模式

```julia
# LP → HiGHS / Gurobi / Clp
function solve_via_highs(problem::OptimizationProblem)
    # 构建 HiGHS 模型
    # 求解
    # 返回统一格式结果
end

# QP → OSQP / Gurobi / Quadprog
function solve_via_osqp(problem::OptimizationProblem)
    # 构建 OSQP 模型
    # 求解
    # 返回统一格式结果
end

# NLP → Ipopt / KNITRO
function solve_via_ipopt(problem::OptimizationProblem)
    # 构建 Ipopt 问题
    # 求解
    # 返回统一格式结果
end
```

---

## 12. 实施路线图

### 阶段 1：核心框架（M1-M2）

| 里程碑 | 交付内容 | 依赖 | 预计工作量 |
|---|---|---|---|
| M1: 问题表示 | `OptimizationProblem` 类型、`OptimizationContract`、`DerivativeContract` | 无 | ~200 行 |
| M2: 问题分析 | `analyze()`、结构检测、问题分类 | M1 | ~300 行 |
| M3: 求解器注册表 | `SolverCapability`、`register_solver!`、`matching_solvers` | M2 | ~200 行 |
| M4: 路由引擎 | `OptimizationRoutePolicy`、`optimize_plan()`、资格门 | M2+M3 | ~400 行 |
| M5: 执行框架 | `solve()` 主循环、回退、验证 | M4 | ~300 行 |
| M6: 证书体系 | `OptimizationCertificate`、`OptimalityVerification` | M5 | ~200 行 |

### 阶段 2：算法实现（M3-M4）

| 里程碑 | 交付内容 | 依赖 | 预计工作量 |
|---|---|---|---|
| M7: 无约束算法 | 最速下降、LBFGS、Newton (CG) | M5 | ~800 行 |
| M8: 约束算法 | SQP、Active Set LP/QP | M5+M7 | ~1200 行 |
| M9: 最小二乘 | LM、Gauss-Newton | M5 | ~500 行 |
| M10: 内点法 | IPM LP → IPM QP → IPM NLP | M5+M6+内层求解器 | ~1000 行 |
| M11: 内层集成 | KKT 系统预条件器、委托测试 | M6+M10 | ~400 行 |

### 阶段 3：外部后端子（M4-M5）

| 里程碑 | 交付内容 | 依赖 |
|---|---|---|
| M12: Ipopt 桥接 | NLP → Ipopt 封装 | M1-M6 |
| M13: OSQP 桥接 | QP → OSQP 封装 | M1-M6 |
| M14: HiGHS/Gurobi 桥接 | LP/MILP → 外部求解器 | M1-M6 |
| M15: AD 集成 | ForwardDiff/Zygote/FiniteDiff 自动检测与适配 | M2 |

### 阶段 4：可观测性与 V&V

| 里程碑 | 交付内容 | 依赖 |
|---|---|---|
| M16: 遥测 | 优化指纹、历史存储、自适应建议 | M6 |
| M17: 测试套件 | 验证矩阵（68+ 问题）、回归测试 | M7-M15 |
| M18: 物理集成 | Vlasov-Fokker-Planck / Maxwell 验证 | M7-M17 |

---

## 13. 与现有项目的集成关系

| 现有项目 | 与本框架的关系 | 集成方式 |
|---|---|---|
| `AdaptiveLinearSolvers.jl` | 内层线性求解器 | 委托调用，继承类型/路由/证书体系 |
| `king-mixture-approximation` | 本框架的早期用户 | King mixture 的 MLE/正则化作为用例 |
| `julia-optimization-ecosystem` | 生态调研与参考 | 提供算法选择和设计参考 |
| `optimization-theory-and-techniques` | 优化理论调研 | 提供理论基础和文献支撑 |
| `Scientific.jl` / `Optimization.jl` | 可选外部桥接 | 通过后端能力模型注册 |

---

## 14. 开放问题

1. **矩阵自由优化的自动微分契约**：当目标仅为函数句柄时，通过自动微分（AD）自动提供梯度，若 AD 不可用则回退到有限差分。但有限差分的精度和成本应如何纳入路由决策？

2. **KKT 预条件器的自动选择**：SQP/IPM 中 KKT 系统的预条件器对收敛至关重要。是否要在自动路由中包含"预条件器测试"步骤（短预算试运行多个预条件器）？

3. **非凸问题的全局性保证**：框架是否应集成多起点、模拟退火或全局搜索作为无梯度方法的回退？当前设计只保证局部最优性。

4. **混合整数问题的 NLP 松弛**：MINLP 的外逼近和 B&B 需要调用 NLP 求解器作为子问题。NLP 求解器的选择是否应进入路由策略，还是由 MINLP 算法内部硬编码？

5. **整数变量与分类**：混合整数问题的指纹和历史自适应策略是否需要与连续问题完全不同的路由逻辑？当前的统一指纹可能不够表达离散结构。