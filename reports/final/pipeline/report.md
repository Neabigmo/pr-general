# PRI-General 单文件流程输出

最终 pre-Boltz 队列：**102,600 条**；2,287 个不同蛋白；73,272 个不同 RNA。

本次构建来源：`raw_source_candidates`。未运行 Boltz，也未修改 data/release 中的原始或阶段文件。

## 输出文件

- `PRI-General_pre_boltz_detailed.csv`：最终队列，每行一条 pair，保留原始字段并附加长度、GC、hash、质量、cluster、degree 和 benchmark 审计字段。
- `general_stage_audit.csv`：源级候选记录的逐阶段进入/未进入标记。
- `source_ingestion_summary.csv`：四类来源实际读取行数及来源级过滤结果。
- `raw_build_vs_release_comparison.csv`：raw-to-final 结果与 release General 的独立对照，不作为输入。
- `stage_waterfall.csv`：从输入到最终队列的数量瀑布。
- `protein_degree.csv`、`rna_degree.csv`：最终队列中每个唯一蛋白/RNA 的配对次数。
- `field_completeness.csv`、`quality_summary.csv`：字段完整性和质量统计。

## 运行方式

`python scripts/pri_general_full_pipeline.py --root H:/2026try/8.15PR/1DATA`

默认从 `data/archive/legacy_pipeline/core/candidates/` 开始重建；若要下载可直接下载的来源快照，在命令后加入 `--download`。MMseqs2 会优先使用项目归档的 `archive/legacy_pipeline/tools/mmseqs/bin/mmseqs.exe`。只有显式加入 `--from-release-snapshot` 才会使用 release General 复核。
