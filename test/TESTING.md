# 测试组织与执行规则

## 默认回归集

`Pkg.test()` 必须运行所有仍受支持功能的测试。当前默认入口为 `test/runtests.jl`，按功能加载以下测试组：

| 目录 | 覆盖范围 |
|---|---|
| `test/direct/` | 直接分解、稀疏输入与残差语义 |
| `test/planning/` | 路线规划、数学资格与路线锁定 |
| `test/krylov/` | GMRES 与 FGMRES 的执行报告 |
| `test/preconditioner/` | 固定/可变预条件器的资格与分流 |
| `test/preconditioner_builders/` | Jacobi 构造、版本缓存与实际 Krylov 使用 |
| `test/operators/` | 矩阵自由与块耦合算子 |
| `test/diagnostics/` | 条件信息版本、诊断预算与数值状态区分 |
| `test/spectral/` | Hermitian 算子的预算化 Lanczos 谱诊断 |
| `test/resources/` | 资源准入、预算终止与后端能力登记 |
| `test/telemetry/` | 指纹、历史、证书与拒绝状态 |

功能目录不是版本目录。新的仍受支持行为必须加入相应功能目录；不得为了保存旧版本测试而创建“已通过测试集合”。历史版本由 Git tag 完整保存。

## 版本复现

Git tag 固定当时的源码、测试和 `Project.toml`。例如复现 `v0.0.2`：

```powershell
git checkout v0.0.2
julia --project=. -e 'using Pkg; Pkg.test()'
```

发布时还必须在 `docs/releases/vX.Y.Z.md` 记录测试命令、Julia 与关键依赖版本、通过数和已知边界。

## 昂贵测试层级

默认测试只包含快速单元与集成测试。昂贵测试按用途放入下列目录，并由环境变量显式选择：

| 层级 | 目录 | PowerShell 执行方式 |
|---|---|---|
| 默认 | 默认功能目录 | `julia --project=. -e 'using Pkg; Pkg.test()'` |
| Nightly | `test/nightly/` | `$env:ALS_TEST_TIER = 'nightly'; julia --project=. test/runtests.jl` |
| Manual | `test/manual/` | `$env:ALS_TEST_TIER = 'manual'; julia --project=. test/runtests.jl` |
| 全部 | 两个昂贵层级 | `$env:ALS_TEST_TIER = 'all'; julia --project=. test/runtests.jl` |

`nightly` 用于大规模 Vlasov--Maxwell 算例、性能基准、GPU 或 MPI 等可自动调度但不宜每次运行的测试。`manual` 用于需要人工确认、受限硬件或外部环境的测试。当前两个目录仅建立入口约定，尚未加入这类昂贵算例。

## 废弃与覆盖率

只有 API 被正式废弃或其语义不再支持时，才可从默认测试移除。移除时必须同时提供变更说明、迁移说明和可定位的 Git 历史。

代码覆盖率只能辅助发现未执行代码，不能替代回归测试。特别是数学资格门、残差语义、GMRES/FGMRES 分流及块耦合行为，必须使用具有明确断言的功能测试验证。
