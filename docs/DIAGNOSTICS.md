# 条件数与数值诊断

本模块只提供诊断和路线排序信息，不改变数学资格门。条件数大不能证明 Hermitian 或正定，历史或诊断结果也不能绕过 CG、MINRES、GMRES 与 FGMRES 的既有资格条件。

## 默认行为

`ConditioningPolicy()` 默认只验证调用方提供的 `ConditioningInfo`。它不计算完整条件数、不额外分解矩阵、不增加算子应用，也不把“未知”标为“病态”。

外部信息必须至少说明指标、作用算子、矩阵版本、来源、证据和可靠性。例如：

```julia
info = ConditioningInfo(
    estimate=2.0e8,
    metric=:kappa_2,
    operator=:original,
    matrix_version=:implicit_step_42,
    source=:assembly_diagnostic,
    reliable=true,
    evidence=:external,
)
problem = AdaptiveLinearProblem(A, b;
    conditioning=info,
    matrix_version=:implicit_step_42,
)
```

版本不匹配返回 `:stale`，不可靠信息返回 `:untrusted`，无效指标或作用算子返回 `:invalid`。只有已接受的、表明近秩亏的外部信息才会使显式矩阵的直接路线优先考虑 `SVD`。

## 受预算估计

估计必须显式启用，并且只适用于预算内的小型显式稠密方阵：

```julia
policy = ConditioningPolicy(
    estimation=:full,
    budget=DiagnosticBudget(max_seconds=1.0, max_matrix_dimension=256),
)
diagnosis = diagnose(problem; policy=policy)
```

`estimation=:cheap` 计算一范数倒数条件数；`estimation=:full` 计算二范数条件数。超出时间或维度预算、稀疏矩阵和矩阵自由算子只返回未执行原因，不进行隐式转换或完整条件数前置计算。

## 区分失败模式

`NumericalDiagnosis` 保存互不混淆的字段：

| 字段 | 示例状态 | 含义 |
|---|---|---|
| `conditioning.state` | `:ill_conditioned`、`:near_rank_deficient` | 已接受条件信息的数值解释 |
| `iteration_state` | `:stagnated`、`:breakdown`、`:nonconverged` | Krylov 后端的可观察状态 |
| `preconditioner_state` | `:suspected_failure` | 预条件器与停滞或 breakdown 相关的提示，不是因果证明 |
| `overall_state` | `:iterative_stagnation`、`:near_rank_deficient` | 面向日志与路线回退的摘要 |

通过 `TelemetryPolicy(level=:basic, output=OutputRequest(diagnostics=true))` 获取按需诊断；设置 `OutputRequest(certificate=true)` 时，诊断也会写入 `RouteCertificate`。
