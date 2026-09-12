# 大规模条件与谱诊断

## Lanczos 谱区间

对已由 `MathematicalContract.hermitian` 认证的方阵，调用方可以显式请求 Lanczos 谱区间估计：

```julia
policy = ConditioningPolicy(
    spectral_estimation=:lanczos,
    spectral_steps=20,
    budget=DiagnosticBudget(
        max_seconds=2.0,
        max_operator_applications=20,
    ),
)
diagnosis = diagnose(problem; policy=policy)
```

该过程只使用 `mul!(y, A, x)`，因此同时适用于显式稀疏矩阵与 `MatrixFreeOperator`。它需要正的时间和算子调用预算；缺失 Hermitian 证据、预算或方阵条件时不会执行，并在 `NumericalDiagnosis.notes` 中记录原因。

`SpectralInfo.lambda_min` 与 `lambda_max` 是 Lanczos Ritz 区间的近似值，不是严格特征值证书。对于已认证正定系统，足够小的估计下端点会标记 `spectral_state=:near_singular`；这是一项数值诊断，不会自行宣布矩阵秩亏或放宽/改变方法资格。

## 与条件数诊断的关系

`ConditioningPolicy.estimation` 的显式稠密条件数估计与 `spectral_estimation=:lanczos` 相互独立：前者面向小型显式方阵，后者面向大型 Hermitian 或矩阵自由问题。默认两者均关闭。

一般非正规矩阵若没有伴随作用和额外求解预算，不存在可靠且通用的低成本完整条件数估计。此时应保留外部 `ConditioningInfo`，并结合迭代残差、停滞、breakdown 和预条件器状态解释求解表现。
