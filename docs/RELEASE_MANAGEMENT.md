# 发布说明管理规则

本规则约束 `AdaptiveLinearSolvers` 的版本说明与 GitHub Release 发布。

## 1. 两层说明

每个已发布版本必须同时具有下列两种说明，二者服务不同读者，不能相互替代。

| 层级 | 位置 | 读者与用途 | 长度 |
|---|---|---|---|
| 简要说明 | `CHANGELOG.md` 对应版本节与 GitHub Release 正文 | 使用者快速判断是否升级 | 3--7 条要点 |
| 详细说明 | `docs/releases/vX.Y.Z.md` | 开发者、数值方法使用者与审计 | 结构化完整记录 |

GitHub Release 正文必须是简要说明，且必须链接到同一标签下的详细说明。不得将详细说明全文直接复制到 GitHub Release 正文。

## 2. 简要说明的固定内容

GitHub Release 与 `CHANGELOG.md` 应以相同的事实范围描述版本。简要说明必须包含：

1. 版本号与发布日期。
2. 一段不超过两句的版本定位。
3. `Added`、`Changed`、`Fixed` 或 `Breaking changes` 中实际发生的类别；无对应事实时省略。
4. 测试结论，包括测试总数、通过数以及 Julia 版本。
5. 指向 `docs/releases/vX.Y.Z.md` 的固定标签链接。

数学资格、路线适用条件、接口参数、已知限制、性能边界、迁移细节和验证环境的完整记录只放入详细说明。

## 3. 详细说明的最低结构

`docs/releases/vX.Y.Z.md` 至少包含：

1. 发布范围与版本定位。
2. 按模块说明新增或修改的行为。
3. 公共接口和兼容性影响。
4. 数学或数值方法的适用边界。
5. 测试命令、环境、通过结果与未覆盖项。
6. 已知限制和下一阶段范围。

详细说明是版本的可追溯技术记录；发布后不得通过静默修改既有标签中的文件来改写其事实。

## 4. 发布任务顺序

1. 完成一个经测试的批次，更新 `Project.toml`、`CHANGELOG.md` 与详细说明。
2. 在本地创建发布提交；版本号仅在批次完成时更新，不随单个微小提交更新。
3. 运行该批次约定的测试，并将结果写入简要与详细说明。
4. 创建带注释标签 `vX.Y.Z`，推送 `main` 与标签。
5. 创建 GitHub Release：正文写入简要说明，并链接该标签中的详细说明。
6. 若发布说明有误，优先修订 GitHub Release 正文并另行提交更正记录；不得重写既有发布标签。

## 5. GitHub Release 正文模板

```markdown
## vX.Y.Z

一句话说明本批版本定位。

### Added

- 核心新增能力。

### Fixed

- 已修复的使用者可见问题。

### Verification

- `Pkg.test()`：N/N 通过，Julia X.Y.Z。

详细的接口、数学边界、测试环境和未覆盖项见：
[vX.Y.Z 详细发布说明](https://github.com/tangwang-USTC/AdaptiveLinearSolvers/blob/vX.Y.Z/docs/releases/vX.Y.Z.md)
```

