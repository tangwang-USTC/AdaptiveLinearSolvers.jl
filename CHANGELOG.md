# Changelog

本文件记录面向使用者的版本变化，格式遵循 Keep a Changelog 的思想。

## [Unreleased]

### Added (Phase 0+1)

### Added (Phase 2)

- 新增 `:bunchkaufman` 直接路线：对称不定矩阵的 Bunch-Kaufman LDLt 分解（`bunchkaufman(Hermitian(A))`），填补 Cholesky（需 SPD）与 LU（通用但无视对称性）之间的数学资格间隙。Hermitian 非正定矩阵可直接求解。
- 新增 `:lsqr` / `:lsmr` 迭代路线：最小二乘 Krylov 求解器（Krylov.jl），支持方形和矩形线性系统，无需 Hermitian/SPD 资格。通过 `Lock(:lsqr)` / `Lock(:lsmr)` 显式请求。
- 为 `CountingOperator` 添加 `adjoint` 支持，使 Krylov.jl 双乘积方法（`A * v` + `A' * u`）可被框架适配。
- 测试覆盖新增 19 项（Bunch-Kaufman 资格 + 求解验证、LSQR/LSMR 方形与矩形求解、矩形稠密 QR 改进测试），总通过数 178 → 197。

### Fixed (Phase 0+1)

- 新增 `:direct` 路线作为 `A\b` 通配后备直接求解器，在 `_DIRECT_CAPABILITIES`、`_direct_order` 和 `_eligibility` 中注册。
- 新增对 IterativeSolvers.jl 后端 `reltol`/`abstol` 关键字参数的适配（替代已废弃的 `tol`）。
- 显著扩展测试覆盖：直接路线中的 QR/SVD/:direct 求解验证（9 项）、稀疏 SPD 和稀疏矩形矩阵测试。
- Krylov 求解验证（CG/MINRES/BiCGStab）及其迭代控制边界条件测试集（4 项边界测试）。
- 单位预条件器构建与求解验收测试。
- 新增 fallback 测试组：回退链成功、资格拒绝、预算终止、SPD 合约路由选择。
- 规划层测试：QR Lock 孤立、组合 Lock/Prefer 策略、可变 PC → FGMRES 替代、family Lock + iterative Lock 空集。
- 资源与边界值测试：自动后端选择、未知后端名、GPU 请求、全 `ResourceBudget` 默认值、无内存预算资源评估。

### Fixed

- 修正 `select_backend` 中 `:stdlib` 后端未注册 `:direct` 路线的问题（返回 `:no_implemented_backend_for_route`）。
- 修正 IterativeSolvers.jl v0.9.4 `tol=` → `reltol=`/`abstol=` 的 API 变更。
- 修正超定矩形方程组 QR 测试：改为构造列空间内真实解验证精确性。
- 修正 fallback 测试中 `notes` 字段访问（默认 telemetry :off 时 `certificate=nothing`）。
- 修正 `Prefer(:gmres)` + 可变预条件器组合策略断言：`planned_routes` 包含 FGMRES 及直接回退路线。
- 修正 CUDA 后端测试断言：执行模式不匹配返回 `:no_implemented_backend_for_route` 而非 `:backend_adapter_not_implemented`。`

## [0.0.4] - 2026-09-08

### Added

- 增加版本匹配的 `ConditioningInfo`、受预算的 `ConditioningPolicy` 与独立 `NumericalDiagnosis`，区分病态、近秩亏、停滞、breakdown、预条件器可疑失效和未收敛。
- 增加字段可控的 `FingerprintProfile`、相似历史检索、路线成功率建议、预条件器复用提示以及显式受信任本机历史持久化。
- 增加 `trace`、`diagnostic` 遥测级别和 `TelemetryBudget`，只采样后端已有残差历史。
- 增加 `ResourceBudget`、`BackendPolicy` 与后端能力登记；标准库直接法和 `Krylov.jl` 在串行 CPU 上可执行，GPU、MPI 与外部适配器显式标记为未实现。
- 将默认回归测试按功能拆分，并增加诊断、历史和资源预算测试。

### Changed

- 默认测试入口继续运行所有支持功能；当前发布批次通过 `85/85` 项测试。

## [0.0.3] - 2026-09-08

### Added

- 接入 `Krylov.jl` 作为无预条件器 CG、MINRES、GMRES、FGMRES 与 BiCGStab 的执行后端，并返回统一的迭代报告和预算终止状态。
- 增加可应用预条件器算子接口，将固定预条件器映射到 Krylov 左预条件器 `M`，将可变预条件器映射到 FGMRES 右预条件器 `N`。
- 增加 `MatrixFreeOperator`，使未显式组装矩阵的线性作用可自动分流到 Krylov 路线。
- 在理论文档中补充矩阵自由算子的定义、Vlasov-Fokker-Planck/Maxwell 例子与路线边界，并由设计文档链接引用。
- 增加 `BlockLayout` 与 `BlockOperator`，保留多字段耦合结构并支持块内显式/矩阵自由算子组合。

### Fixed

- 为自定义矩阵自由与块算子补齐 Krylov 所需的 `size(A, dimension)` 协议。

## [0.0.2] - 2026-09-08

### Added

- 增加无数值执行副作用的 `plan(problem, policy)` 规划层、路线能力和分层判定记录。
- 将直接路线执行器改为只消费 `RoutePlan.execution_routes`，并扩展规划、稀疏、零右端项和秩亏测试。
- 增加迭代法与预条件器数学资格接口，明确 GMRES 与 FGMRES 对可变预条件器的分流规则。
- 将迭代资格接入规划结果，区分数学合格但尚无后端的 `unavailable_routes` 与实际可执行路线。

## [0.0.1] - 2026-09-08

### Added

- 建立数学自适应线性求解器的理论基线。
- 建立分层路线覆盖设计，支持 `Auto`、`Prefer`、`Lock` 和 `Forbid`。
- 建立代码、文档、测试与发布共用版本的 Git/GitHub 治理规则。
- 建立 `0.0.1` Julia 数学核心：直接路线、数学资格门、分层路线控制与默认关闭的基础/指纹遥测。
- 包含 `SolveStatus`、残差验收、按需 `RouteCertificate`、本地开发安装说明和 `Pkg.test()` 测试目标。

### Changed

- 明确 FGMRES（Flexible Generalized Minimum Residual，柔性广义最小残量法）仅在同一次线性求解内预条件器可变、非线性或其固定性未知时替代 GMRES。
- 将条件数定义为带范数、算子、证据和矩阵版本的可选输入；补充按规模的估计与未知条件信息路线。
- 定义默认关闭的 `off`、定制字段 `basic`、标签驱动 `fingerprint`、`trace` 与 `diagnostic` 五级遥测，以及由标签指纹和事件触发、受数学资格约束的历史建议机制。
