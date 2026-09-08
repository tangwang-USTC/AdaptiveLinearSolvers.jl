# Changelog

本文件记录面向使用者的版本变化，格式遵循 Keep a Changelog 的思想。

## [Unreleased]

## [0.1.3] - 2026-09-08

### Changed

- 将公共 `solve(A, b; kwargs...)` 入口前置到路由源文件，核心实现保持在其后。

## [0.1.2] - 2026-09-08

### Added

- 增加 `SolveStatus`、绝对/相对残差验收策略与按需 `RouteCertificate`。
- 对资格拒绝和全部直接路线数值失败返回结构化状态，而非仅抛出异常。

## [0.1.1] - 2026-09-08

### Fixed

- 声明 Julia 标准库 `Test` 为测试目标，使 `Pkg.test()` 能构建测试环境。
- 说明未注册本地包的 `Pkg.develop(path=...)` 加载方式与包级测试命令。

## [0.1.0] - 2026-09-08

### Added

- 建立数学自适应线性求解器的理论基线。
- 建立分层路线覆盖设计，支持 `Auto`、`Prefer`、`Lock` 和 `Forbid`。
- 建立代码、文档、测试与发布共用版本的 Git/GitHub 治理规则。
- 建立 `0.1.0` Julia 数学核心：直接路线、数学资格门、分层路线控制与默认关闭的基础/指纹遥测。

### Changed

- 明确 FGMRES（Flexible Generalized Minimum Residual，柔性广义最小残量法）仅在同一次线性求解内预条件器可变、非线性或其固定性未知时替代 GMRES。
- 将条件数定义为带范数、算子、证据和矩阵版本的可选输入；补充按规模的估计与未知条件信息路线。
- 定义默认关闭的 `off`、定制字段 `basic`、标签驱动 `fingerprint`、`trace` 与 `diagnostic` 五级遥测，以及由标签指纹和事件触发、受数学资格约束的历史建议机制。
