# AdaptiveLinearSolvers.jl

一个面向 Julia 的、以数学资格为边界的自适应线性代数求解路由器。

项目的工作名称为 `AdaptiveLinearSolvers.jl`。在注册 Julia General registry 或创建 GitHub 远程仓库前，名称、UUID、许可证和组织归属仍可调整。

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

当前仓库仅建立理论和设计基线，尚未实现 Julia 源程序，尚未添加外部求解器依赖，也尚未配置 GitHub 远程。
