# 审计表占位说明

`general_stage_audit.csv` 当前约 2.26 GB，按项目上传策略不提交到 GitHub。

它是当前发布流程的逐行审计表，不是运行代码所需的源输入。完整文件应保留在本地发布归档的同一路径：

```text
data/final/audit/general_stage_audit.csv
```

缺少该文件时，依赖逐行审计的来源阶段总览和部分分析图不能从仓库副本重建；Parquet 阶段快照和最终索引仍可按仓库中实际存在的文件使用。补回该文件后再运行 `run_general.py validate`，验证入口才会完整通过。
