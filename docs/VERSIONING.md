# 版本、Git 与远程治理

## 术语与缩略词

| 缩略词 | 首次定义 |
|---|---|
| API | Application Programming Interface，应用程序编程接口 |
| HTTPS | Hypertext Transfer Protocol Secure，安全超文本传输协议 |

## 1. 版本对象

仓库中的理论文档、设计文档、源程序、测试、基准和说明文档属于同一个发布对象。每个 Git tag 固定一个完整、可重建的项目快照；文档不单独漂移为无法对应代码的“最新说明”。

发布版本采用 `MAJOR.PHASE.BATCH`，即“主版本.阶段.小批量改动版”：

| 版本变化 | 适用情形 |
|---|---|
| `MAJOR` | 公共 API、数学契约语义、证书格式或缓存语义发生不兼容变化 |
| `PHASE` | 在同一主版本内完成一个兼容能力阶段，例如新增一组路线、后端扩展、预条件器、证据类型或文档章节 |
| `BATCH` | 已计划并验收的兼容改动批次，例如一组实现、文档、测试或数值判定修复；单个提交、纯排版或代码排布调整不单独递增 |

在 `0.x.x` 阶段，`PHASE` 可以包含试验性 API 的不兼容调整；从 `1.0.0` 起严格遵循稳定 API 承诺。

Git commit 是细粒度开发与审计单位；版本号是可复现发布批次单位。只有当一组有明确目标的改动已完成、文档同步并通过该批次规定的测试后，才递增 `BATCH`。纯重构、排版和单点修复只创建提交，并并入下一次批量版本。

`Project.toml` 中的包版本、`CHANGELOG.md` 的发布标题和 Git tag 必须一致，例如三者均为 `0.1.0` / `v0.1.0`。

## 2. 文档版本

- `README.md`：当前项目入口与范围。
- `THEORY.md`：数学假设、资格规则和可验证结论。
- `DESIGN.md`：软件架构、API 语义与实现边界。
- `docs/`：用户说明、迁移指南、发布记录和维护规则。
- `CHANGELOG.md`：面向用户的版本变化，不替代 Git 提交历史。

文档与代码在同一分支、同一提交中演进。发布后，Git tag 是该版本文档的权威快照；未来若建立 GitHub Pages，应发布 `stable` 指向最新 tag、`dev` 指向 `main`，但不复制维护第二套文档。

## 3. Git 分支与提交

远端默认主分支为 `main`，只接收已审查、可复现的集成提交。项目开发分支采用：

```text
agent/p/adaptive-linear-solvers/<topic>
```

示例：

```text
agent/p/adaptive-linear-solvers/router-core
agent/p/adaptive-linear-solvers/krylov-extension
agent/p/adaptive-linear-solvers/docs-api
agent/p/adaptive-linear-solvers/petsc-extension
```

提交信息使用 Conventional Commits 的简洁子集：

```text
feat(router): add eligibility gate for Hermitian methods
fix(cache): invalidate numeric factorization on version change
docs(theory): define evidence levels for positive definiteness
test(krylov): cover variable-preconditioner fallback
```

一个拉取请求只处理一种可描述的变更。理论规则改变必须同时修改相应设计、测试和变更日志；实现改变若影响公开行为，必须同时修改设计或用户说明。

## 4. GitHub 远程

该仓库作为独立 GitHub 项目管理。远程仅使用 HTTPS：

```text
origin  https://github.com/<owner>/AdaptiveLinearSolvers.jl.git
```

首次公开推送前应确认：

- GitHub 所有者或组织。
- 仓库名称是否采用 `AdaptiveLinearSolvers.jl`。
- 私有或公开可见性。
- 开源许可证；在未确定前不应声明可再分发。
- 是否启用 issue、discussions、GitHub Actions 和 GitHub Pages。

建议在 GitHub 上保护 `main`：要求拉取请求、至少一次审查、通过测试、禁止强推。发布由受保护 `main` 创建带注释 tag `vMAJOR.PHASE.BATCH`，再生成 GitHub Release；未经明确授权不得自动推送、创建 Release 或修改远程设置。

## 5. 发布最小清单

1. `Project.toml`、`CHANGELOG.md` 和文档中的版本号一致。
2. 理论契约、设计、源程序、测试和示例在同一提交中一致。
3. 运行测试与基准，保存命令、Julia 版本、依赖版本和结果摘要。
4. 提交到 `main` 后创建 `vMAJOR.PHASE.BATCH` tag。
5. 创建 GitHub Release，并链接对应 tag 的文档快照。

## 6. 初始仓库状态

当前版本为 `0.0.3`：主版本 `0` 表示试验期，阶段 `0` 表示初始基线尚未进入独立能力阶段，小批量改动版 `3` 归并 Krylov 执行后端、可应用预条件器、矩阵自由算子和块耦合算子。`0.0.3` 已在 Julia `1.12.6`、Krylov `0.10.9` 下通过 `48/48` 包测试，耗时约 `15.7 s`。后续零散修改仅记录为 Git commit；满足批量发布阈值后再递增版本。GitHub 远程、许可证与公开发布仍需在项目身份确认后单独执行。
