# AdaptiveLinearSolvers.jl 差距分析与后续任务规划

## 版本与状态

| 项目 | 值 |
|---|---|
| 当前基线 | v0.0.4 |
| 审查日期 | 2026-09-13 |
| 源码规模 | 1,725 行 / 14 文件 |
| 测试规模 | 380 行 / 11 组 / 85 项通过 |
| 文档覆盖 | DESIGN.md (251行) + THEORY.md (341行) + ARCHITECTURE.md (500行) + docs/ 6 篇 |

---

## 一、差距分析总表

颜色标记：🔴 未实现 / 🟡 部分实现 / ✅ 已完成

### A. 求解器路线

| 路线 | 资格门 | 执行器 | 测试 | 状态 |
|---|---|---|---|---|
| LU 直接法 | ✅ | ✅ | ✅ | **已实现** |
| Cholesky 直接法 | ✅ | ✅ | ✅ | **已实现** |
| QR 直接法 | ✅ | ✅ | 🔴 无独立测试 | 🟡 功能完整，测试不足 |
| SVD 直接法 | ✅ | ✅ | 🔴 无独立测试 | 🟡 功能完整，测试不足 |
| `A\b` (`:direct`) 通用路线 | 🔴 未注册到 _DIRECT_CAPABILITIES | ✅ 执行器存在 | 🔴 不可达 | 🔴 **死代码** |
| CG 迭代法 | ✅ | ✅ (Krylov.jl) | 🔴 仅资格门测试，无求解验证 | 🟡 功能完整，测试不足 |
| MINRES 迭代法 | ✅ | ✅ (Krylov.jl) | 🔴 仅资格门测试，无求解验证 | 🟡 功能完整，测试不足 |
| GMRES 迭代法 | ✅ | ✅ (Krylov.jl) | ✅ 基本求解验证 | **已实现** |
| FGMRES 迭代法 | ✅ | ✅ (Krylov.jl) | ✅ 基本求解验证 | **已实现** |
| BiCGStab 迭代法 | ✅ | ✅ (Krylov.jl) | 🔴 **完全无测试** | 🟡 功能完整，测试缺失 |
| LSQR/LSMR 最小二乘 | 🔴 未注册 | 🔴 未实现 | 🔴 无 | **未实现** |
| Bunch-Kaufman/LDL | 🔴 未注册 | 🔴 未实现 | 🔴 无 | **未实现** |

### B. 预条件器

| 预条件器 | 构造 | 缓存 | 测试 | 状态 |
|---|---|---|---|---|
| Jacobi | ✅ | ✅ | ✅ | **已实现** |
| Identity | ✅ | ✅ | 🔴 无独立测试 | 🟡 功能完整，测试不足 |
| ILU (不完全 LU) | 🔴 未实现 | 🔴 | 🔴 | **未实现**（DESIGN.md §14 提及） |
| AMG (代数多重网格) | 🔴 未实现 | 🔴 | 🔴 | **未实现**（DESIGN.md §14 提及） |
| 块 Schur 补 | 🔴 未实现 | 🔴 | 🔴 | **未实现**（DESIGN.md §16 提及） |
| 物理预条件器 | 🔴 未实现 | 🔴 | 🔴 | **未实现**（DESIGN.md §16 提及） |

### C. 后端适配器

| 后端 | 登记 | 实现 | 测试 | 状态 |
|---|---|---|---|---|
| stdlib (LinearAlgebra) | ✅ | ✅ 所有直接法 | ✅ | **已实现** |
| krylov (Krylov.jl) | ✅ | ✅ 全部 5 种方法 | 🟡 部分 | **已实现** |
| iterativesolvers | ✅ | 🟡 仅 GMRES | ✅ GMRES | 🟡 部分实现 |
| linearsolve | ✅ | 🟡 仅 GMRES，无预条件器 | ✅ GMRES | 🟡 部分实现 |
| petsc | ✅ | 🔴 未实现 | 🔴 | **未实现** |
| cuda | ✅ | 🔴 未实现 | 🔴 | **未实现** |
| mpi | ✅ | 🔴 未实现 | 🔴 | **未实现** |

### D. 矩阵自由与块算子

| 功能 | 实现 | 测试 | 状态 |
|---|---|---|---|
| MatrixFreeOperator (正向作用) | ✅ | ✅ | **已实现** |
| BlockOperator (块耦合) | ✅ | ✅ | **已实现** |
| CountingOperator (计数+预算) | ✅ | ✅ | **已实现** |
| 伴随作用 A† | 🔴 未实现 | 🔴 | **未实现**（THEORY.md §3.1） |
| 块算子预条件器 (块三角/Schur) | 🔴 未实现 | 🔴 | **未实现**（DESIGN.md §16） |
| Jv/Jtv 产品接口 | 🔴 未实现 | 🔴 | **未实现** |

### E. 数值诊断

| 诊断项 | 实现 | 测试 | 状态 |
|---|---|---|---|
| 外部条件信息验证 | ✅ | ✅ | **已实现** |
| Lanczos 谱估计 | ✅ | ✅ | **已实现** |
| 迭代状态诊断 (停滞/breakdown) | ✅ | ✅ | **已实现** |
| 预条件器状态诊断 | ✅ | ✅ | **已实现** |
| 条件数估计 (cheap: rcond_1) | ✅ | ✅ | **已实现** |
| 条件数估计 (full: kappa_2) | ✅ | ✅ | **已实现** |
| 一范数逆估计器 | 🔴 未实现 | 🔴 | **未实现**（THEORY.md §3.2） |
| 后向误差计算 | 🔴 未实现 | 🔴 | **未实现**（DESIGN.md §5） |

### F. 文档

| 文档 | 状态 |
|---|---|
| DESIGN.md | ✅ **已有**（251行） |
| THEORY.md | ✅ **已有**（341行） |
| ARCHITECTURE.md | ✅ **已有**（500行，v0.1 2026-09-12） |
| docs/BACKENDS.md | ✅ **已有** |
| docs/PRECONDITIONERS.md | ✅ **已有** |
| docs/SPECTRAL_DIAGNOSTICS.md | ✅ **已有** |
| docs/OBSERVABILITY.md | ✅ **已有** |
| docs/RELEASE_MANAGEMENT.md | ✅ **已有** |
| CHANGELOG.md | ✅ **已有** |

### G. 测试覆盖缺口

| 缺口 | 影响 | 当前覆盖 |
|---|---|---|
| CG 求解验证测试 | 未验证 CG 实际收敛 | 仅资格门 |
| MINRES 求解验证测试 | 未验证 MINRES 实际收敛 | 仅资格门 |
| BiCGStab 完全无测试 | 未验证 BiCGStab 资格 + 求解 | **零覆盖** |
| QR/SVD 路线求解测试 | 未验证直接法完整覆盖 | 只有::cholesky 和 ::lu 在测试中 |
| 多路线回退链测试 | 未验证连续失败→回退→BudgetTerminated | 仅有单一回退场景 |
| 混合资格+策略集成测试 | 未验证多种 contract+policy 组合 | 基本 |
| 稀疏矩阵扩展测试 | 稀疏矩阵只测了 `sparse(A)` 基本 | 1 行 |
| 大规模/性能测试 | nightly/manual 目录为空 | **零覆盖** |
| 物理集成 (M5) 测试 | 无 Vlasov-Fokker-Planck/Maxwell 参考 | **零覆盖** |

---

## 二、任务规划

### 阶段 0 — 代码修复与补漏（Pre-v0.1.0）

此阶段修复已有代码中的缺陷，不增加新能力。

| 编号 | 任务 | 文件 | 类型 |
|---|---|---|---|
| H-01 | 将 `:direct` 路线注册到 `_DIRECT_CAPABILITIES`，使 Julia `\` 通用路线可通过 plan() 访问 | `planning.jl` L52-57 | 修复 |
| H-02 | 为 `:qr`、`:svd` 补充独立求解测试，确保残差验收正确 | `test/direct/runtests.jl` | 测试 |
| H-03 | 为 `:identity` 预条件器构造补充独立测试 | `test/preconditioner_builders/runtests.jl` | 测试 |
| H-04 | 修复 `backends.jl` 中缺失的 `_validate_iteration_control` 边界条件测试（`max_iterations=0` 和 `max_seconds=0`） | `test/krylov/runtests.jl` | 测试 |

### 阶段 1 — 测试覆盖补齐（Priority: High）

补齐已有功能的测试覆盖，不做代码修改。

| 编号 | 任务 | 预计新增测试 | 测试文件 |
|---|---|---|---|
| T-01 | **BiCGStab 资格门 + 求解验证测试**（Hermitian 不定矩阵，验证收敛） | ~15 行 | `test/krylov/runtests.jl` |
| T-02 | **CG 求解验证测试**（SPD 矩阵，带预条件器/不带预条件器） | ~20 行 | `test/krylov/runtests.jl` |
| T-03 | **MINRES 求解验证测试**（Hermitian 不定矩阵） | ~15 行 | `test/krylov/runtests.jl` |
| T-04 | **回退链集成测试**（多路线连续失败 → FallbackSuccess / BudgetTerminated） | ~30 行 | `test/planning/runtests.jl` 新建 `test/fallback/runtests.jl` |
| T-05 | **混合资格+策略组合测试**（部分 Lock + 部分 Auto + 强/弱 contract 组合） | ~25 行 | `test/planning/runtests.jl` |
| T-06 | **稀疏矩阵扩展测试**（多种稀疏模式，使用 SparseArrays） | ~20 行 | `test/direct/runtests.jl` |
| T-07 | **资源预算边界值测试**（max_seconds=0, max_iterations=0, max_memory_bytes=0 等异常输入） | ~15 行 | `test/resources/runtests.jl` |

**阶段 1 通过标准**：新增 ~140 行测试，回归后总通过数从 85 项增至 ~115 项。

### 阶段 2 — 求解器路线扩展（Priority: High）

| 编号 | 任务 | 预计工作量 | 依赖 |
|---|---|---|---|
| R-01 | **添加 `:bunchkaufman` / `:ldl` 路线**：对称不定直接分解（扩展 `_solve_route` + `_eligibility`） | ~30 行源码 + ~20 行测试 | 阶段 0 |
| R-02 | **添加 LSQR/LSMR 路线**：最小二乘 Krylov 求解器（扩展 `_ITERATIVE_CAPABILITIES` + `_solve_route`；Krylov.jl 已支持） | ~35 行源码 + ~25 行测试 | 阶段 0 |

### 阶段 3 — 预条件器扩展（Priority: Medium）

| 编号 | 任务 | 预计工作量 | 依赖 |
|---|---|---|---|
| P-01 | **ILU 预条件器构造器**：稀疏矩阵的 Incomplete LU 分解（新构造器 + 资格测试） | ~100 行源码 + ~30 行测试 | 阶段 1 |
| P-02 | **不完全 Cholesky 预条件器**：SPD 稀疏矩阵（可选，依赖 ILU 框架） | ~60 行源码 + ~20 行测试 | P-01 |
| P-03 | **AMG 预条件器**：AlgebraicMultigrid.jl 桥接（扩展能力模型） | ~80 行源码 + ~30 行测试 | 阶段 1 |

### 阶段 4 — 数值诊断与矩阵自由增强（Priority: Medium）

| 编号 | 任务 | 预计工作量 | 依赖 |
|---|---|---|---|
| D-01 | **后向误差计算**：`η = ‖b - Ax‖/(‖A‖‖x‖+‖b‖)` 加入 `_residual_metrics` 和 RouteCertificate | ~30 行源码 + ~15 行测试 | 阶段 1 |
| D-02 | **一范数逆估计器**：利用 `A\e` 估计 ‖A⁻¹‖₁（受 DiagnosticBudget 约束） | ~50 行源码 + ~20 行测试 | 阶段 1 |
| D-03 | **MatrixFreeOperator 伴随接口**：`apply_adjoint!` 可选字段，用于 LSQR/LSMR 和条件估计 | ~40 行源码 + ~20 行测试 | R-02 |
| D-04 | **块算子基本预条件器**：块对角 Jacobi（利用 BlockLayout 的分块信息） | ~60 行源码 + ~25 行测试 | 阶段 2 |

### 阶段 5 — 后端增强（Priority: Low-Medium）

| 编号 | 任务 | 预计工作量 | 依赖 |
|---|---|---|---|
| B-01 | **IterativeSolvers 扩展**：增加 CG、MINRES 适配（当前只有 GMRES） | ~60 行源码 + ~20 行测试 | 阶段 1 |
| B-02 | **LinearSolve 增强**：支持预条件器映射、多种算法选择 | ~80 行源码 + ~25 行测试 | 阶段 1 |
| B-03 | **nightly 测试第一个用例**：中等规模 (~1000×1000) 稀疏矩阵基准 | ~40 行测试 | 阶段 1 |

### 阶段 6 — Vlasov-Fokker-Planck/Maxwell 集成验证（M5）（Priority: Low）

| 编号 | 任务 | 预计工作量 | 依赖 |
|---|---|---|---|
| V-01 | **Vlasov 矩方程线性系统代表性测试**（小规模合成矩阵，块结构 + 不定 + 矩阵自由） | ~60 行测试 + 理论说明 | 阶段 2+4 |
| V-02 | **Maxwell curl-curl 系统代表性测试**（SPD 半定 + 矩阵自由） | ~50 行测试 + 理论说明 | 阶段 2+4 |

---

## 三、阶段测试策略

每个阶段遵循统一验证流程：

```
阶段实施
  ├── 实现阶段功能（源码变更）
  ├── 写功能测试（验证新行为的正确性）
  ├── 回归运行：Pkg.test() 确认 0 失败
  ├── 边界测试：负面输入、错误条件、空值/零值
  └── 更新 CHANGELOG.md（记录变更）
```

### 测试类型定义

| 测试类型 | 定位 | 覆盖要求 |
|---|---|---|
| 单元测试 | 单函数/单方法，纯数学 | 每个资格门、每个诊断状态、每个路线 |
| 集成测试 | solve() 端到端 | 每个路线至少一个成功案例 |
| 边界测试 | 异常输入、零值、空值 | `max_*=0`、`nothing`、`PropertyEvidence()` |
| 回退测试 | 多路线失败链 | 主线→回退→BudgetTerminated |
| 回归测试 | `Pkg.test()` 默认运行 | 阶段实施结束时 0 失败 |
| 昂贵测试 | `ALS_TEST_TIER=nightly` | 大规模、长时间、性能基准 |

### 阶段验收条件

| 阶段 | 代码增量 | 测试增量 | 验收命令 |
|---|---|---|---|
| 阶段 0 | ≤10 行 | ~50 行 | `julia --project=. -e 'using Pkg; Pkg.test()'` 通过 ≥ 100 项 |
| 阶段 1 | 0 行 | ~140 行 | 通过 ≥ 115 项 |
| 阶段 2 | ~65 行 | ~45 行 | 通过 ≥ 125 项 |
| 阶段 3 | ~240 行 | ~80 行 | 通过 ≥ 140 项 |
| 阶段 4 | ~180 行 | ~80 行 | 通过 ≥ 155 项 |
| 阶段 5 | ~180 行 | ~85 行 | 通过 ≥ 170 项 |
| 阶段 6 | 测试为主 | ~110 行 | 通过 ≥ 180 项 |

---

## 四、工作标准（继承项目规范）

1. **所有变更必须通过 `Pkg.test()` 回归**，禁止提交未通过回归的代码。
2. **数学资格门优先级高于性能启发式**。ILU 即使快但不收敛，也需记录不收敛原因。
3. **新增路线必须填写所有四个组件**：资格门 (`_eligibility` / `qualify_iterative`) + 执行器 (`_solve_route`/`_execute_route`) + 能力登记 (`_DIRECT_CAPABILITIES` / `_ITERATIVE_CAPABILITIES`) + 测试。
4. **后端适配器不能改变核心数学资格门**。IterativeSolvers 即使声称收敛，也要通过本项目残差验收。
5. **CHANGELOG.md 必须同步更新**。格式遵循 Keep a Changelog。
6. **每次发布必须递增版本号**（`Project.toml` 中 `version`），并在 `docs/releases/` 记录测试环境和通过数。
7. **不暴露未授权功能**：不联网、不调用外部 API，除非获得明确授权。

---

## 五、执行顺序建议

```
阶段 0 ──→ 阶段 1 ──→ 阶段 2 ──→ 阶段 3 ──→ 阶段 4 ──→ 阶段 5 ──→ 阶段 6
 (修复)    (测试补齐)  (路线扩展)  (预条件器)  (诊断增强)  (后端)    (物理集成)
                      │                              ↑
                      └── 依赖阶段 1 的测试覆盖 ──────┘
```

**建议优先执行**：阶段 0 → 阶段 1 → 阶段 2 或 阶段 3（可并行）。

预计总计新增代码：~700 行源码 + ~590 行测试，目标总通过数从 85 → 180+ 项。