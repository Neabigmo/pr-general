# PRI-General 独立工作区

这是当前最终版的 PRI-General 普通数据集工作区。这里不包含旧版本、历史报告、旧 PPT 或临时聚类数据库。

筛选边界和当前发布数量集中记录在 `config/pipeline.yaml`；运行入口和分析脚本都从这里读取规则，避免同一参数在多个文件中漂移。

## GitHub 上传范围

`data/final/audit/general_stage_audit.csv` 约 2.26 GB，按上传策略不提交，仅以 `data/final/audit/README.md` 占位说明。`data/source_candidates/npinter_raw.parquet` 和 `data/final/complexes/PRI-General_final_complexes.csv` 通过 Git LFS 保存；其余当前工作区文件按仓库内容上传。

## 当前数据

数据集按“复合物结构”统计：

- 最终复合物：**101,336 个**；
- 需要预测结构的复合物：**99,576 个**；
- 已有实验结构的 PDB 复合物：**1,760 个**；
- 来源：NPInter 主表 32,079 个、NPInter binding-site 67,185 个、NPInter miRNA 312 个、RCSB PDB 1,760 个；
- 蛋白长度：40–2,000 aa；RNA 长度：10–500 nt；
- Boltz 前预测队列：**92,350 个**。这是筛选队列记录数，和 101,336 个最终复合物结构单位不是同一个统计口径。

## 使用

安装依赖：

```bash
python -m venv .venv
.venv\\Scripts\\activate
python -m pip install -r requirements.txt
```

检查当前最终数据：

```bash
python run_general.py validate
```

重新生成最终复合物索引、来源阶段总览和全部当前统计图：

```bash
python run_general.py analyze
```

从源候选表重新执行筛选流程：

```bash
python run_general.py rebuild
```

`rebuild` 默认读取 `data/source_candidates/`。如果要重新计算 Protein90、Protein30 和 RNA90 聚类，需要在命令中额外提供目标环境的 MMseqs2 参数。

## 目录

```text
PRI-General_final/
├─ run_general.py                    # 独立入口
├─ config/pipeline.yaml             # 当前筛选参数和数量
├─ data/source_candidates/           # 四个来源的当前候选表
├─ data/reference/                  # 普通源快照和 benchmark 参考集
├─ data/final/complexes/             # 一行一个最终复合物的结构索引
├─ data/final/pre_boltz/             # D1、D2、D3 和预测队列
├─ data/final/audit/                 # 当前最终流程的逐行审计表
├─ reports/final/                    # 当前最终统计与审计结果
├─ deliverables/analysis/            # PNG、PDF 和关键 CSV 图表
├─ scripts/                          # 构建、验证和分析脚本；common.py 放共享规则/身份 helper
└─ docs/                             # 简明数据说明
```

PDB 的统计单位是唯一 `pdb_id`；NPInter 来源的统计单位是唯一 `sequence_hash`。底层审计表中的技术行不直接作为复合物数量。
