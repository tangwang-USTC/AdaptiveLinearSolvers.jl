# Changelog

本文件记录面向使用者的版本变化，格式遵循 Keep a Changelog 的思想。

## [Unreleased]

### Added

- 建立数学自适应线性求解器的理论基线。
- 建立分层路线覆盖设计，支持 `Auto`、`Prefer`、`Lock` 和 `Forbid`。
- 建立代码、文档、测试与发布共用版本的 Git/GitHub 治理规则。

### Changed

- 明确 FGMRES（Flexible Generalized Minimum Residual，柔性广义最小残量法）仅在同一次线性求解内预条件器可变、非线性或其固定性未知时替代 GMRES。
- 将条件数定义为带范数、算子、证据和矩阵版本的可选输入；补充按规模的估计与未知条件信息路线。
- 定义默认关闭的 `off`、`summary`、`trace` 与 `diagnostic` 四级遥测，以及由标签指纹和事件触发、受数学资格约束的历史建议机制。
