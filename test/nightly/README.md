# Nightly 测试

此目录只放置耗时较长、需要较大计算资源或依赖专门环境的测试，例如中等规模 Vlasov--Maxwell 算例、性能基准、GPU 和 MPI 验证。

根测试入口只在 `ALS_TEST_TIER=nightly` 或 `ALS_TEST_TIER=all` 时加载本目录直接包含的 `.jl` 文件。每个文件应自行建立 `@testset`，不得依赖执行顺序。
