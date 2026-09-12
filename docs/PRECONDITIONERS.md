# 预条件器构造与版本缓存

## 当前可执行构造器

`build_preconditioner(problem)` 当前只实现以下构造器：

| `PreconditionerBuildPolicy.kind` | 适用输入 | 产生的逆作用 |
|---|---|---|
| `:jacobi` | 显式方阵且对角元不接近零 | `Diagonal(inv.(diag(A)))` |
| `:identity` | 任意方阵 | 显式恒等逆作用 |

返回的 `PreconditionerBuildReport.contract` 是可直接交给 `with_preconditioner` 的 `PreconditionerContract`。随后现有 GMRES、FGMRES、CG、MINRES 与 BiCGStab 路径继续按照固定性、线性性、Hermitian 性和正定性资格门选择方法。

```julia
cache = PreconditionerCache(8)
report = build_preconditioner(problem; cache=cache)
prepared = with_preconditioner(problem, report.contract)
solution = solve(prepared; policy=RoutePolicy(iterative=Lock(:gmres)))
```

## 缓存边界

缓存键包含 `label`、`matrix_version`、矩阵尺寸和构造策略。没有 `matrix_version` 时，构造器仍可建立预条件器，但不会读取或写入缓存。调用方必须在离散算子、网格、时间步装配或物理参数改变后更新 `matrix_version`。

缓存只复用已构造的逆作用，不自动判断物理模型是否仍可复用，也不跨进程持久化。`clear_preconditioner_cache!` 只清理缓存引用，不修改调用方已有的预条件器对象。

## 后续构造器

ILU（Incomplete LU，不完全 LU 分解）、不完全 Cholesky、AMG（Algebraic Multigrid，代数多重网格）、块 Schur 补和物理预条件器尚未实现。它们必须在独立适配批次中明确其矩阵格式、并行能力、参数、失败语义和 GMRES/FGMRES 固定性契约。
