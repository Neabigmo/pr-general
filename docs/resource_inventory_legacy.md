# PRI-General：已有资源、快速补齐项与当前缺口

## 一、包内已经具备的内容

| 内容 | 包内位置 | 当前用途 |
|---|---|---|
| NPInter v5 主表 | `data/source_candidates/npinter_raw.parquet` | 普通 protein–RNA 相互作用候选 |
| NPInter v5 binding-site | `data/source_candidates/npinter_bindingsite_raw.parquet` | RNA 结合位点窗口候选 |
| NPInter v5 miRNA | `data/source_candidates/mirna_protein_raw.parquet` | miRNA–protein 候选 |
| RCSB PDB 候选表 | `data/source_candidates/pdb_raw.parquet` | PDB 来源记录和结构候选信息 |
| 最终复合物索引 | `data/final/complexes/PRI-General_final_complexes.parquet` | 101,336 个最终复合物结构单位 |
| PDB 最终索引 | `data/final/complexes/pdb_final_complexes.csv` | 1,760 个 PDB 实验结构复合物 |
| D1/D2/D3 阶段 | `data/final/pre_boltz/` | 当前筛选阶段快照 |
| 蛋白聚类 | `data/final/pre_boltz/clusters/protein_cluster_assignments.parquet` | 已有 Protein90/Protein30 相关分配 |
| RNA 聚类 | `data/final/pre_boltz/clusters/rna_cluster_assignments.parquet` | 已有 RNA90 相关分配 |
| 统计和代码 | `scripts/`、`reports/final/` | 当前结果的复核和再分析 |

## 二、可以较快补上的内容

### 1. PDB 结构文件

从 `data/final/complexes/pdb_final_complexes.csv` 读取 1,760 个唯一 `pdb_id`，逐个下载 RCSB 的 biological assembly mmCIF，并保留下载日志。随后解析：

1. assembly 中的 protein/RNA 链；
2. 分辨率和实验方法；
3. 原子坐标；
4. 蛋白主链 N/Cα/C/O 完整性；
5. RNA 糖–磷酸骨架原子完整性；
6. 蛋白–RNA 界面残基和核苷酸。

这一步是最重要的补充，因为当前包只有 PDB 候选表和部分历史解析结果，没有完整的当前 mmCIF 坐标包。

### 2. UniProt、RNAcentral 和 Rfam

- 用 UniProt accession 补齐蛋白序列、物种和蛋白名称；
- 用 RNAcentral/miRBase accession 补齐 RNA 序列；
- 用 Rfam family/clan 对 RNA 做家族注释。

原项目归档中已有部分缓存和 Rfam 文件，可作为恢复入口，但需要记录数据库版本、下载日期和映射成功率，不能直接把旧缓存当成完整当前快照。

目前核对到的归档入口包括：

- UniProt：`H:\\2026try\\8.15PR\\1DATA\\data\\archive\\legacy_pipeline\\core\\raw\\uniprot_cache.tsv`；
- RNAcentral/miRNA：`H:\\2026try\\8.15PR\\1DATA\\data\\archive\\legacy_pipeline\\core\\raw\\mirbase\\mirbase_rnacentral.fasta.gz`；
- Rfam：`H:\\2026try\\8.15PR\\1DATA\\data\\archive\\legacy_pipeline\\core\\raw\\rfam\\`；
- RNAInter：`H:\\2026try\\8.15PR\\1DATA\\data\\archive\\legacy_pipeline\\core\\raw\\rnainter\\`；
- IntAct：`H:\\2026try\\8.15PR\\1DATA\\data\\archive\\legacy_pipeline\\core\\raw\\intact\\`；
- ENCODE/RBNS：`H:\\2026try\\8.15PR\\1DATA\\data\\archive\\legacy_pipeline\\expansion\\source_material\\raw\\encode_rbns\\`。

### 3. RNAInter、IntAct、ENCODE/RBNS

原项目归档中有 RNAInter 压缩包、IntAct 压缩包以及 ENCODE/RBNS 元数据或文件清单。快速补齐方式是：下载或恢复原始文件，统一成标准字段，再单独保存为外部证据或 assay 资产。它们不能不加区分地直接并入天然复合物主表。

### 4. Protein40、RNA80 和更完整的家族注释

当前有 Protein90、Protein30 和 RNA90 的已有分配。若需要 Protein40、RNA80，应在目标环境安装 MMseqs2 后重新聚类，并记录阈值、覆盖率、版本和运行命令。当前环境中没有可直接调用的 `mmseqs` 命令。

### 5. Boltz-2 pilot

当前机器可确认的 GPU 是 NVIDIA GeForce RTX 4070 Laptop GPU，约 8 GB 显存。快速补齐方式是：

1. 安装并固定 Boltz-2 代码版本或 commit；
2. 下载并记录模型权重位置；
3. 固定 MSA 工具和数据库版本；
4. 先运行 1,000 条 pilot；
5. 记录显存、失败原因、耗时和输出文件位置。

本包当前没有 Boltz-2 命令、权重或 MSA 数据库，因此尚未具备直接预测条件。

## 三、当前根本没有的内容

以下内容在当前工作区和已核对的原项目材料中没有找到可直接使用的完整版本：

- RNApedia 数据快照；
- POSTAR3 数据快照；
- RNA3DB 结构组件映射；
- 准确对应的 PRDB v3 数据包；
- 完整 RNAcentral 快照；
- 已固定版本的 Boltz-2 权重和 MSA 数据库；
- 已对当前 1,760 个 PDB 复合物完成的统一 interface backbone 原子完整性表。

其中，原项目归档里的 `rnaprodb`、`rnaprotdb` 等相关材料不等同于 PRDB v3，不能替代它。

原项目中还存在 PDB 的历史结构解析目录
`H:\\2026try\\8.15PR\\1DATA\\data\\archive\\legacy_pipeline\\core\\experimental_structures\\`，其中包括结构元数据、接触对和部分 assembly 解析结果；但这些文件没有随本包复制，也没有证明它们覆盖当前全部 1,760 个 PDB 复合物。

## 四、使用时需要注意

当前包可以支持 PRI-General 的数据查看、数量核对、筛选阶段复核和已有统计重算，但不能宣称：

- 所有 1,760 个 PDB 复合物的坐标都已在包内；
- 所有来源都已经纳入 PRI-General；
- Rfam、RNA3DB 或其他 RNA family 注释已经完成；
- Boltz-2 已经安装并可以直接运行。

因此，这个包代表的是“当前已经具备且已整理出的内容”，不是包含所有计划数据源的完整二期资源包。
