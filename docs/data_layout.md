# 当前文件说明

## 最终数据

| 文件 | 数量 | 含义 |
|---|---:|---|
| `data/final/complexes/PRI-General_final_complexes.parquet` | 101,336 | 一行一个最终复合物，PDB 保留一个 PDB 的全部不同组件 |
| `data/final/complexes/pdb_final_complexes.csv` | 1,760 | 最终保留的 PDB 实验结构复合物 |
| `data/final/pre_boltz/prediction_queue.parquet` | 92,350 | 当前 Boltz 前预测队列 |
| `data/final/pre_boltz/general_d1_hard_qc.parquet` | 150,093 | D1 后队列 |
| `data/final/pre_boltz/general_d2_exact_dedup.parquet` | 150,093 | exact 去重阶段快照 |
| `data/final/pre_boltz/general_d3_local_rna_representatives.parquet` | 143,414 | 同蛋白内 RNA90 代表样本 |

## 来源和参考数据

`data/source_candidates/` 保存四个来源的当前候选表，用于从头重跑。`data/reference/benchmark/` 只用于检查训练泄漏，不属于 PRI-General 最终数据。

`data/final/audit/general_stage_audit.csv` 是当前流程的逐行审计表，保留每条来源记录在各阶段是否通过。它很大，但属于当前最终流程的追溯数据。

GitHub 副本不提交这个约 2.26 GB 的 CSV；原目录中的 `data/final/audit/README.md` 是占位说明。需要完整审计追溯时，应从本地发布归档补回该文件。

## 统计结果

统计表和图位于 `reports/final/analysis/` 与 `deliverables/analysis/`。主统计单位是复合物结构数，不是底层技术记录数。
