# PRI-General 是什么

PRI-General 收集 protein–RNA 复合物候选，用于后续结构预测和序列建模。当前版本把所有来源统一成一套 protein–RNA 序列记录，再按复合物结构进行统计。

数据主要来自：

- NPInter 主表：普通的 protein–RNA 相互作用；
- NPInter binding-site：RNA 结合位点窗口；
- NPInter miRNA：miRNA 相关相互作用；
- RCSB PDB：已有实验结构的 protein–RNA 复合物。

处理顺序是：读取来源候选表 → 解析蛋白和 RNA 序列 → 检查字符和长度 → 合并重复序列组合 → 进行 Protein90、Protein30 和 RNA90 聚类 → 对同一蛋白的相近 RNA 做代表性保留 → 执行蛋白和家族限额 → 得到最终队列。

当前最终结果包含 101,336 个复合物结构单位，其中 1,760 个来自 PDB、已有实验结构；其余 99,576 个是后续需要预测结构的复合物候选。
