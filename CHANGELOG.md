# Changelog

本文件记录面向使用者的版本变化，格式遵循 Keep a Changelog 的思想。

## [Unreleased]

## [0.0.2] - 2026-09-08

### Added

- 增加无数值执行副作用的 `plan(problem, policy)` 规划层、路线能力和分层判定记录。
- 将直接路线执行器改为只消费 `RoutePlan.execution_routes`，并扩展规划、稀疏、零右端项和秩亏测试。
- 增加迭代法与预条件器数学资格接口，明确 GMRES 与 FGMRES 对可变预条件器的分流规则。
- 将迭代资格接入规划结果，区分数学合格但尚无后端的 `unavailable_routes` 与实际可执行路线。
- 接入 `Krylov.jl` 作为无预条件器 CG、MINRES、GMRES、FGMRES 与 BiCGStab 的执行后端，并返回统一的迭代报告和预算终止状态。
- 增加可应用预条件器算子接口，将固定预条件器映射到 Krylov 左预条件器 `M`，将可变预条件器映射到 FGMRES 右预条件器 `N`。
- 增加 `MatrixFreeOperator`，使未显式组装矩阵的线性作用可自动分流到 Krylov 路线。

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
