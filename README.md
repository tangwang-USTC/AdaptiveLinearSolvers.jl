# AdaptiveLinearSolvers.jl

一个面向 Julia 的、以数学资格为边界的自适应线性代数求解路由器。

项目的工作名称为 `AdaptiveLinearSolvers.jl`。在注册 Julia General registry 或创建 GitHub 远程仓库前，名称、UUID、许可证和组织归属仍可调整。

## 术语与缩略词

| 缩略词 | 首次定义 |
|---|---|
| CG | Conjugate Gradient，共轭梯度法 |
| MINRES | Minimum Residual，最小残量法 |
| AMG | Algebraic Multigrid，代数多重网格 |
| ILU | Incomplete LU，不完全 LU 分解 |
| PETSc | Portable, Extensible Toolkit for Scientific Computation，可移植可扩展科学计算工具包 |
| GPU | Graphics Processing Unit，图形处理器 |

## 项目目标

- 为稠密、稀疏、结构化和矩阵自由线性系统提供统一接口。
- 依据可审计的数学证据选择直接法、Krylov 方法和预条件器。
- 支持部分强制：锁定某一决策层，其他层保持自动选择。
- 返回解、数值验证结果和完整的路由证书。
- 通过 Julia package extensions 接入 `Krylov.jl`、`LinearSolve.jl`、AMG、ILU、PETSc 和 GPU 后端。

## 非目标

- 不在核心层编码 Vlasov-Fokker-Planck、矩方程或 Maxwell 的物理知识。
- 不把未经认证的“近似对称”或“经验上正定”当作 CG、Cholesky 或 MINRES 的许可。
- 不以无约束地试遍所有求解器作为自适应策略。

## 文档

- [理论基础](THEORY.md)
- [软件设计](DESIGN.md)
- [版本、Git 与远程治理](docs/VERSIONING.md)
- [变更日志](CHANGELOG.md)

## 计划中的目录

```text
src/                 Julia 源程序
test/                纯数学测试与回归基准
docs/                使用说明、版本与发布记录
examples/            最小可复现示例
benchmarks/          受版本控制的数学基准定义
```

## 当前状态

`0.1.0` 已提供无外部依赖的数学核心：显式稠密/稀疏矩阵的直接路线、数学契约资格门、分层路线控制、可选条件数证据，以及默认关闭的 `off`、`basic`、`fingerprint` 遥测。

这一阶段刻意未实现 Krylov 迭代法、预条件器、矩阵无显式表示、GPU、MPI 或外部后端；未实现的路线会明确拒绝，而不会伪装为可用能力。

包入口为 `src/AdaptiveLinearSolvers.jl`。测试已写入 `test/runtests.jl`；本次仅建立代码，尚未执行测试。
