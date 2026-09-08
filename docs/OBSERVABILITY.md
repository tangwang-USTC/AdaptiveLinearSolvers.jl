# 可观测性与历史自适应

## 默认边界

`TelemetryPolicy(level=:off)` 仍是默认值：不构造指纹、不写入历史、不持久化，也不为遥测额外计算真实残差、条件数或谱信息。`FingerprintProfile` 逐字段控制 `representation`、`size_band`、`structure`、`conditioning` 与 `execution`；关闭的字段以 `nothing` 返回，且不参与相似性匹配。

## 历史建议

`HistoryStore` 是容量受限的内存记录库。`similar_records` 只匹配当前指纹中启用的字段与相同的 `label`；`route_advice` 统计候选路线的历史成功率，并给出预条件器复用提示。

当 `plan` 或 `solve` 接收 `history` 时，建议只可重排已经通过资格门的自动路线。任何 `Lock` 或 `Prefer` 选择保持原有顺序；历史记录不能使未认证的 SPD 系统使用 CG，也不能使可变预条件器绕过 FGMRES。

预条件器复用只返回 `:reuse_candidate`、`:no_reuse_evidence` 或 `:not_applicable`，不自动复用对象。实际缓存键、生命周期和物理模型版本仍由调用方控制。

## 受控持久化

`save_history(path, store)` 与 `load_history!(store, path)` 只在调用方明确指定路径时读写历史。格式使用 Julia `Serialization`，因此只应加载本机受信任路径中的文件；它不是跨 Julia 版本或跨语言的长期数据交换格式。

## Trace 与 diagnostic

`TelemetryPolicy(level=:trace)` 和 `:diagnostic` 只保留后端已产生的残差历史，不会自行重新计算残差。`TelemetryBudget(max_trace_samples=N)` 限制保留样本数，`sample_every` 控制抽样步长；零样本预算只报告是否截断。

`level=:diagnostic` 允许同时输出诊断对象，但条件数估计仍必须由 `ConditioningPolicy` 显式开启并受其独立预算约束。遥测级别不能隐式提升诊断开销。
