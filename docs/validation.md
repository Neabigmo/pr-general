# 当前版本核查

运行（校验值从 `config/pipeline.yaml` 读取）：

```bash
python run_general.py validate
```

当前核查项目：

- 最终复合物索引为 101,336 行；
- PDB 复合物为 1,760 个；
- 预测候选为 99,576 个；
- Boltz 前队列为 92,350 条；
- D1、D2、D3 阶段分别为 150,093、150,093、143,414 条；
- `complex_id`、预测队列 `sample_id` 和 `sequence_hash` 唯一；
- 预测队列蛋白长度、RNA 长度和 benchmark 冲突检查通过；
- `sequence_hash` 与蛋白、RNA 序列重新计算结果一致。

这里的 1,760 个 PDB 是筛选后保留的唯一 `pdb_id` 数量。原始 PDB 候选表有 4,993 个唯一 PDB 条目，二者不是同一个数字。PDB 坐标文件没有被虚构为已经下载；当前工作区保存的是 PDB 来源记录和最终复合物索引。
