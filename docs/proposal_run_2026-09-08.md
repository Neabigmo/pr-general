# 四主来源方案试跑记录（2026-09-08）

本次试跑把来源处理和统一筛选拆开：只有同时具备 `protein_sequence`、`rna_sequence`、`exact_protein_mapping=true`、`exact_rna_mapping=true` 的记录才进入 exact pair 池；同一 pair 的多条证据保留在 evidence 表中，但 pair 只计一次。

## 可复现流程

```powershell
F:\anaconda3\python.exe scripts/download_encode_eclip.py --output-dir data/downloads/encode_eclip --workers 8
F:\anaconda3\python.exe scripts/external_resolver_yield.py --root . --external-source-root <local-source-cache> --sample-rows 100000 --output-dir reports/final/rebuild/external_resolver_yield
F:\anaconda3\python.exe scripts/run_proposal.py --root . --external-source-root <local-source-cache> --output-dir reports/final/rebuild/proposal_run
```

`<local-source-cache>` 应包含 RNAInter、IntAct 和历史 ENCODE/RBNS 目录。运行输出目录位于 `reports/final/rebuild/`，该目录被 `.gitignore` 排除；本文件和 eCLIP 元数据清单是可提交的审计摘要。

## NPInter / PDB 漏斗

| 输入表 | raw | protein exact | RNA exact | 双侧 exact | 通过字母表/长度 | unique pair |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| NPInter main | 564,646 | 557,036 | 338,585 | 334,378 | 41,244 | 41,080 |
| NPInter binding-site | 197,052 | 196,748 | 197,052 | 196,748 | 155,739 | 114,963 |
| NPInter miRNA | 316 | 313 | 316 | 313 | 312 | 312 |
| RCSB PDB | 41,759 | 41,759 | 41,759 | 41,759 | 17,400 | 5,392 |

NPInter main 的主要修正是沿用现有流水线的 RNA `T→U` 标准化。修正后，main 不再是 0 条可用；剩余主要损失来自 RNA 长度范围，而不是双侧 resolver 全部失败。

## 外部 resolver 抽样

每个归档抽取 100,000 行，仅使用 accession/xref，不使用基因名猜序列。

| 来源 | 过滤结果 | 实验关系 | 蛋白序列解析 | RNA 序列解析 | 双侧 exact |
| --- | ---: | ---: | ---: | ---: | ---: |
| RNAInter RH | 0 RNA–protein | 0 | 0 | 0 | 0 |
| RNAInter RP | 88,129 RNA–protein | 28,315 | 0 | 0 | 0 |
| IntAct | 147 protein–RNA | 52 direct/non-physical | 2 | 0 | 0 |

这是相对于当前项目本地 ID→序列表的保守下界；RNAInter/IntAct 仍需补充 RefSeq、Ensembl、RNAcentral 等专用 resolver 后再判断最终 yield。

## ENCODE eCLIP

原有 RBNS/RNAcompete 矩阵不作为 interaction pair 输入。重新下载官方已发布的 eCLIP `bed narrowPeak` 双生物学重复 `[1,2]` 文件：

- 477 个 peak BED 文件；GRCh38 252 个，hg19 225 个。
- 压缩体积约 40.2 MB；逐文件 MD5 校验通过。
- 当前只完成 peak 资产下载，尚未完成基因组坐标、转录本、链方向到 RNA 窗口，以及 RBP 蛋白唯一序列映射，因此暂不进入 pair 池。

## 当前结果

- 证据记录：214,695
- exact sequence pair：161,747
- 唯一蛋白：2,946
- 唯一 RNA：122,759
- 30% 单来源限制：改为 `OBSERVE`，暂不作为硬门禁。
- Protein40/RNA80/Rfam 尚未生成，严格家族隔离未完成。
- Boltz-2 可执行文件、权重和 checkpoint 不在本机，因此本次止于预测前候选池。

官方 eCLIP 资产检索使用 ENCODE REST API；eCLIP 流程文档说明双重复可复现峰以 merged peak BED 形式提供：

- https://www.encodeproject.org/help/rest-api/
- https://www.encodeproject.org/documents/1f171ac6-a36a-41ac-b632-741aeb47aad2/@@download/attachment/eCLIP_analysisSOP_v2.3.pdf
