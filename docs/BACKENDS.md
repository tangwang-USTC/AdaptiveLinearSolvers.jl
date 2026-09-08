# 计算资源与后端

## 统一资源预算

`ResourceBudget` 是每次 `solve` 的显式资源上限：

| 字段 | 含义 |
|---|---|
| `max_seconds` | 与 `IterationControl.max_seconds` 取较小值 |
| `max_iterations` | 与 Krylov 迭代上限取较小值；零表示不追加限制 |
| `max_memory_bytes` | 对当前后端已知的最小工作区作保守准入；未知时拒绝有界内存请求 |
| `max_operator_applications` | 当前 Krylov 后端以迭代上限作保守约束，不宣称精确计数 |
| `execution` | `:serial_cpu`、`:gpu` 或 `:distributed_cpu` |

默认预算不额外限制时间、迭代或内存，并且只允许 `:serial_cpu`。请求无法验证的内存上限、未实现的执行形态或超出下限估计的工作区时，求解返回 `BudgetTerminated` 并在路线证书中记录原因。

## 后端能力登记

| 后端标识 | Julia 包 | 当前状态 | 执行形态 |
|---|---|---|---|
| `:stdlib` | `LinearAlgebra` | 已实现 | 串行 CPU 直接法 |
| `:krylov` | `Krylov.jl` | 已实现 | 串行 CPU Krylov 法 |
| `:iterativesolvers` | `IterativeSolvers.jl` | 已登记，未实现 | 串行 CPU |
| `:linearsolve` | `LinearSolve.jl` | 已登记，未实现 | 串行 CPU、GPU |
| `:petsc` | `PETSc.jl` | 已登记，未实现 | 串行 CPU、分布式 CPU |
| `:cuda` | `CUDA.jl` | 已登记，未实现 | GPU |
| `:mpi` | `MPI.jl` | 已登记，未实现 | 分布式 CPU |

`backend_capabilities()` 返回完整登记表。`select_backend` 仅选择“已实现、支持该路线并匹配执行形态”的后端；`BackendPolicy(requested=:petsc)` 等请求会得到 `:backend_adapter_not_implemented`，不会回退到未声明的后端。

## 外部适配器接入原则

接入 IterativeSolvers.jl、LinearSolve.jl、PETSc.jl、CUDA.jl 或 MPI.jl 前，必须分别完成：

1. 后端依赖、版本兼容与许可证审查。
2. `AbstractLinearOperator`、显式稀疏矩阵或分布式矩阵的输入适配。
3. 残差、停止状态、迭代数、时间和内存信息向统一结果对象的映射。
4. 预条件器语义与 GMRES/FGMRES 资格规则的交叉测试。
5. 对应设备或分布式环境的独立回归与性能测试。

后端适配器不能改变本项目的数学资格门、残差验收语义或预算终止含义。
