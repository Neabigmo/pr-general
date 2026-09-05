#!/usr/bin/env python
"""PRI-General 单文件构建、审计与完整 CSV 导出脚本。

这个文件把当前项目中 PRI-General 的可复现流程集中到一个入口：

1. 记录并可选下载来源快照；
2. 默认读取下载/解析后的源级候选表，而不是读取最终 General 表；
3. 执行 D1 硬质量检查、D2 精确 pair 去重、D3 同蛋白内 RNA90 去冗余；
4. 执行单蛋白配对限额和 Protein30 家族限额；
5. 输出按当前输入和规则实时计算的 pre-Boltz 队列完整 CSV；
6. 输出逐阶段审计表、字段完整性、来源、protein/RNA degree 和质量统计。

默认模式使用 `data/source_candidates/` 中的源级候选表，
并重新做全流程。`data/reference/general_source_snapshot.parquet` 只作为参考结果用于对照，
不是默认构建输入。只有显式加入 `--from-release-snapshot` 时，才会读取已经整理好的
General 表；这个模式用于复核，不代表从源数据重建。

源级候选表本身是下载快照经过各来源解析器得到的中间产物：它们保留了无法解析序列、
非法长度和其他失败记录，便于统计每一步损失。源文件下载与解析规则在本文件的
`SOURCE_DOWNLOAD_PLAN` 和 `load_raw_source_candidates` 中明确记录。
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import subprocess
import tempfile
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from common import RAW_SOURCE_DATABASES, load_pipeline_rules, pair_hash


SCRIPT_VERSION = "1.2.0-complex-structure-build"
AA_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")
RNA_RE = re.compile(r"^[ACGU]+$")
AMBIGUOUS_AA = set("XBZJUO")
AMBIGUOUS_RNA = set("RYMKSWHBVDN")

SOURCE_DOWNLOAD_PLAN: list[dict[str, Any]] = [
    {
        "source_id": "npinter_main",
        "source": "NPInter v5",
        "url": "http://bigdata.ibp.ac.cn/npinter5/download/file/interaction_NPInterv5.txt.gz",
        "relative_path": "npinter/interaction_NPInterv5.txt.gz",
        "condition": "只保留 ncRNA-protein binding；按实验方法/数据源分为 E3、E2、E1；纯计算预测 E0 不进入可用 pair。",
    },
    {
        "source_id": "npinter_bindingsite",
        "source": "NPInter v5 binding site",
        "url": "http://bigdata.ibp.ac.cn/npinter5/download/file/bindingsite_NPInterv5.txt.gz",
        "relative_path": "npinter/bindingsite_NPInterv5.txt.gz",
        "condition": "只保留合法基因组区间；转录本/基因可解析、与 exon 有交集；窗口不超过 200 nt。",
    },
    {
        "source_id": "noncode_human",
        "source": "NONCODE v6 human",
        "url": "http://www.noncode.org/datadownload/",
        "relative_path": "noncode/NONCODEv6_human.fa.gz",
        "condition": "作为 NPInter RNA 序列解析来源；仅接受可解析的 NONCODE transcript/gene。具体下载文件由 NONCODE 页面快照确认。",
    },
    {
        "source_id": "noncode_mouse",
        "source": "NONCODE v6 mouse",
        "url": "http://www.noncode.org/datadownload/",
        "relative_path": "noncode/NONCODEv6_mouse.fa.gz",
        "condition": "作为 NPInter RNA 序列解析来源；仅接受可解析的 NONCODE transcript/gene。具体下载文件由 NONCODE 页面快照确认。",
    },
    {
        "source_id": "rnacentral_mirbase",
        "source": "RNAcentral miRBase mirror",
        "url": "https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/sequences/by-database/mirbase.fasta.gz",
        "relative_path": "mirbase/mirbase_rnacentral.fasta.gz",
        "condition": "只保留 NPInter 中 ncType=miRNA 且 RNAcentral/miRBase 与 UniProt 两端都能解析的记录。",
    },
    {
        "source_id": "rcsb_pdb",
        "source": "RCSB PDB",
        "url": "https://search.rcsb.org/rcsbsearch/v2/query",
        "relative_path": "pdb/assembly-1-mmcif/",
        "condition": "实验测定、同时含 protein 和 RNA polymer；排除 ribosome/spliceosome 等大复合体；assembly-1 中 RNA 原子与蛋白原子距离 <=5 Å 且至少 5 个 RNA 原子。",
    },
    {
        "source_id": "uniprot",
        "source": "UniProt REST",
        "url": "https://rest.uniprot.org/uniprotkb/stream",
        "relative_path": "uniprot/uniprot_cache.tsv",
        "condition": "按 NPInter target accession 批量解析 FASTA；缓存 accession→protein sequence；解析失败的 protein 不进入可用 pair。",
    },
]

SYSTEM_MAP = {
    "rbp_general": "sequence_specific_RBP",
    "puf_family": "sequence_specific_RBP",
    "trna_protein": "tRNA_protein",
    "viral_protein_rna": "viral_RNA_protein",
    "rna_processing": "RNA_processing",
    "rna_helicase": "RNA_processing",
    "rna_editing": "RNA_editing",
    "rna_modification_enzyme": "RNA_modification",
    "mirna_machinery": "miRNA_machinery",
    "aptamer_protein": "aptamer",
}
RARE_SYSTEMS = {
    "tRNA_protein", "dsRNA_shape_specific", "bacterial_regulatory_RNP",
    "viral_RNA_protein", "RNA_processing", "RNA_editing", "RNA_modification",
    "miRNA_machinery", "aptamer", "sequence_specific_RBP",
}

RAW_SOURCE_FILES = {
    "npinter_main": "npinter_raw.parquet",
    "npinter_bindingsite": "npinter_bindingsite_raw.parquet",
    "mirna_protein": "mirna_protein_raw.parquet",
    "rcsb_pdb": "pdb_raw.parquet",
}
def project_root(value: Path | None) -> Path:
    return (value or Path(__file__).resolve().parents[1]).resolve()


def nonempty(value: Any) -> bool:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return False
    return str(value).strip().lower() not in {"", "nan", "none", "null"}


def text_seq(value: Any) -> str:
    if not nonempty(value):
        return ""
    return re.sub(r"[^A-Za-z]", "", str(value)).upper()


def normalize_protein(value: Any) -> str:
    return text_seq(value)


def normalize_rna(value: Any) -> str:
    return text_seq(value).replace("T", "U")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_sources(root: Path, download_dir: Path | None = None, force: bool = False) -> list[dict[str, Any]]:
    """按 SOURCE_DOWNLOAD_PLAN 下载可直接下载的文件，并记录所有结果。

    NONCODE 页面、RCSB API 和 UniProt REST 是入口而不是单个静态文件；这些项目
    会被记录为 endpoint_only，不会把网页误保存成 FASTA。实际项目构建使用已有
    归档快照或由对应 parser 解析这些入口。
    """
    destination = (download_dir or root / "data" / "downloads").resolve()
    destination.mkdir(parents=True, exist_ok=True)
    results: list[dict[str, Any]] = []
    endpoint_only = {"noncode_human", "noncode_mouse", "rcsb_pdb", "uniprot"}
    for item in SOURCE_DOWNLOAD_PLAN:
        target = destination / item["relative_path"]
        record = dict(item)
        record["target"] = str(target)
        record["downloaded_at_utc"] = datetime.now(timezone.utc).isoformat()
        if item["source_id"] in endpoint_only:
            record["status"] = "endpoint_only"
            record["note"] = "该来源需要页面选择、API 查询或 accession 批量请求，未把入口网页当作数据文件下载。"
            results.append(record)
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists() and target.stat().st_size > 0 and not force:
            record["status"] = "cached"
            record["bytes"] = target.stat().st_size
            record["sha256"] = sha256_file(target)
            results.append(record)
            continue
        partial = target.with_suffix(target.suffix + ".part")
        try:
            request = urllib.request.Request(item["url"], headers={"User-Agent": "PRI100K-full-pipeline/1.0"})
            with urllib.request.urlopen(request, timeout=300) as response, partial.open("wb") as handle:
                while True:
                    chunk = response.read(1024 * 1024)
                    if not chunk:
                        break
                    handle.write(chunk)
            partial.replace(target)
            record["status"] = "downloaded"
            record["bytes"] = target.stat().st_size
            record["sha256"] = sha256_file(target)
        except Exception as error:  # noqa: BLE001
            if partial.exists():
                partial.unlink()
            record["status"] = "failed"
            record["error"] = repr(error)
        results.append(record)
    manifest = root / "reports" / "final" / "rebuild" / "download_manifest.json"
    manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_text(json.dumps(results, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return results


def _stable_source_prefix(source: str) -> str:
    return {
        "NPInter_v5": "NPI",
        "NPInter_v5_bindingsite": "NBS",
        "NPInter_v5_miRNA": "NPM",
        "RCSB_PDB": "PDB",
    }.get(source, "SRC")


def _raw_source_columns() -> list[str]:
    return [
        "protein_sequence", "rna_sequence", "protein_accession", "rna_accession",
        "protein_name", "rna_name", "organism", "taxonomy_id", "system_class",
        "system_subtype", "rna_class", "source_database", "source_record_id",
        "source_url_or_identifier", "interaction_evidence",
        "experimental_structure_available", "pdb_id", "window_start", "window_end",
        "window_strand", "notes",
    ]


def load_raw_source_candidates(root: Path, candidates_dir: Path, output_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """从源级候选表构建 General 的原始候选池。

    这些 parquet 不是 release/general.parquet 的副本，而是四个来源解析器在
    序列解析之后保留下来的候选记录。函数不在这里删除序列为空的记录；它们
    会进入 D1 并在阶段审计中标为失败。这样“源级候选数”和“质量筛选后数量”
    是可区分、可复核的。
    """
    rules = load_pipeline_rules(root)
    candidates_dir = candidates_dir.resolve()
    frames: list[pd.DataFrame] = []
    audit: list[dict[str, Any]] = []
    print(f"[raw] 读取源级候选目录: {candidates_dir}", flush=True)
    for source_id, filename in RAW_SOURCE_FILES.items():
        path = candidates_dir / filename
        if not path.exists():
            raise FileNotFoundError(f"缺少源级候选文件: {path}")
        expected = RAW_SOURCE_DATABASES[source_id]
        print(f"[raw]   {filename}", flush=True)
        frame = pd.read_parquet(path, columns=_raw_source_columns())
        frame = frame.reset_index(drop=True)
        original_rows = len(frame)
        frame["source_input_row"] = np.arange(1, len(frame) + 1, dtype=np.int64)
        observed = frame["source_database"].fillna("").astype(str)
        source_name_ok = observed.eq(expected)
        if source_id == "npinter_main":
            evidence = frame["interaction_evidence"].fillna("").astype(str)
            source_pass = evidence.str.match(r"^E[1-3](?::|$)")
            source_reason = np.where(source_pass, "", "not_usable_evidence_or_missing_evidence")
        elif source_id == "npinter_bindingsite":
            starts = pd.to_numeric(frame["window_start"], errors="coerce")
            ends = pd.to_numeric(frame["window_end"], errors="coerce")
            window_ok = starts.notna() & ends.notna() & (ends >= starts) & ((ends - starts + 1) <= rules["binding_site_window_max_nt"])
            source_pass = window_ok
            source_reason = np.where(source_pass, "", "invalid_or_overlong_binding_site_window")
        elif source_id == "mirna_protein":
            source_pass = frame["rna_class"].fillna("").astype(str).str.lower().eq("mirna")
            source_reason = np.where(source_pass, "", "not_mirna_record")
        else:
            exp = frame["experimental_structure_available"].fillna(False).astype(bool)
            source_pass = exp & frame["pdb_id"].map(nonempty)
            source_reason = np.where(source_pass, "", "not_an_experimental_pdb_pair")
        source_pass = source_pass & source_name_ok
        source_reason = np.where(~source_name_ok, "unexpected_source_database", source_reason)
        frame["source_filter_pass"] = source_pass.astype(bool)
        frame["source_filter_reason"] = pd.Series(source_reason, index=frame.index).astype(str)
        frame["protein_sequence"] = frame["protein_sequence"].map(normalize_protein)
        frame["rna_sequence"] = frame["rna_sequence"].map(normalize_rna)
        resolved = frame["source_filter_pass"] & frame["protein_sequence"].ne("") & frame["rna_sequence"].ne("")
        unresolved_rows = int((frame["source_filter_pass"] & ~resolved).sum())
        # 不能解析成两端序列的记录仍在 source_ingestion_summary 中统计，
        # 但不带入后续 pair 筛选，避免把空序列当成一个普通 QC 记录并造成无谓内存占用。
        frame = frame.loc[resolved].copy()
        frame["source_filter_pass"] = True
        frame = frame.drop(columns=["source_filter_reason"], errors="ignore")
        prefix = _stable_source_prefix(expected)
        pair_suffix = [pair_hash(p, r)[:12] if p and r else "no_sequence" for p, r in zip(frame["protein_sequence"], frame["rna_sequence"])]
        # source_input_row 让同一来源中重复的 source_record_id 也有唯一主键；
        # D2 后仍按蛋白-RNA 序列 hash 合并真正重复 pair。
        frame["sample_id"] = [
            f"{prefix}_{str(record_id).strip() or 'no_record'}_{suffix}_{row:07d}"
            for record_id, suffix, row in zip(frame["source_record_id"], pair_suffix, frame["source_input_row"])
        ]
        frame["sequence_hash"] = [pair_hash(p, r) for p, r in zip(frame["protein_sequence"], frame["rna_sequence"])]
        frame["benchmark_conflict"] = False
        frame["raw_source_id"] = source_id
        frames.append(frame)
        audit.append({
            "source_id": source_id,
            "source_database": expected,
            "input_file": str(path.relative_to(root)) if path.is_relative_to(root) else str(path),
            "rows_read": int(original_rows),
            "source_filter_pass_rows": int(source_pass.sum()),
            "source_filter_fail_rows": int((~source_pass).sum()),
            "resolved_pair_rows_retained": int(len(frame)),
            "unresolved_sequence_rows": unresolved_rows,
            "nonempty_protein_rows": int(frame["protein_sequence"].ne("").sum()),
            "nonempty_rna_rows": int(frame["rna_sequence"].ne("").sum()),
            "unique_sequence_hash_rows": int(frame["sequence_hash"].nunique()),
        })
    raw = pd.concat(frames, ignore_index=True, sort=False)
    if raw["sample_id"].duplicated().any():
        raise ValueError("源级候选拼接后 sample_id 不唯一；请检查来源行号和来源前缀。")
    audit_frame = pd.DataFrame(audit)
    output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(audit_frame, output_dir / "source_ingestion_summary.csv")
    return raw, audit_frame


def compare_with_release_general(root: Path, raw: pd.DataFrame, final: pd.DataFrame, output_dir: Path) -> pd.DataFrame:
    """把 raw-to-final 结果和当前 release General 做对照，不把 release 当输入。"""
    rules = load_pipeline_rules(root)
    reference_path = root / "data" / "reference" / "general_source_snapshot.parquet"
    rows: list[dict[str, Any]] = []
    if not reference_path.exists():
        rows.append({"metric": "reference_general_exists", "raw_build_value": True, "release_reference_value": False, "status": "无法对照：参考文件不存在"})
        out = pd.DataFrame(rows)
        write_csv(out, output_dir / "raw_build_vs_release_comparison.csv")
        return out
    reference = pd.read_parquet(reference_path, columns=["protein_sequence", "rna_sequence", "source_database"])
    ref_p = reference["protein_sequence"].map(normalize_protein)
    ref_r = reference["rna_sequence"].map(normalize_rna)
    ref_hash = {pair_hash(p, r) for p, r in zip(ref_p, ref_r)}
    raw_hash = set(raw.loc[raw["source_filter_pass"], "sequence_hash"])
    final_hash = set(final["recomputed_sequence_hash_full"])
    rows.extend([
        {"metric": "raw_source_candidate_rows", "raw_build_value": int(len(raw)), "release_reference_value": int(len(reference)), "status": "已验证/参考对照"},
        {"metric": "raw_source_pass_rows", "raw_build_value": int(raw["source_filter_pass"].sum()), "release_reference_value": None, "status": "已验证"},
        {"metric": "raw_source_strict_main_rows", "raw_build_value": int((raw["source_filter_pass"].fillna(False).astype(bool) & raw["protein_sequence"].str.fullmatch(AA_RE).fillna(False) & raw["rna_sequence"].str.fullmatch(RNA_RE).fillna(False) & raw["protein_sequence"].str.len().between(rules["protein_length_min"], rules["protein_length_max"]) & raw["rna_sequence"].str.len().between(rules["rna_length_min"], rules["rna_length_max"])).sum()), "release_reference_value": None, "status": "已验证"},
        {"metric": "raw_source_unique_pairs", "raw_build_value": int(len(raw_hash)), "release_reference_value": int(len(ref_hash)), "status": "已验证/参考对照"},
        {"metric": "reconstructed_final_rows", "raw_build_value": int(len(final)), "release_reference_value": None, "status": "已由当前源数据和规则重新构建；不使用历史队列数量作目标"},
        {"metric": "reconstructed_final_pairs_in_release_general", "raw_build_value": int(len(final_hash & ref_hash)), "release_reference_value": int(len(final_hash)), "status": "已验证"},
        {"metric": "release_general_pairs_found_in_raw_source_pool", "raw_build_value": int(len(ref_hash & raw_hash)), "release_reference_value": int(len(ref_hash)), "status": "已验证"},
        {"metric": "release_general_pairs_not_found_in_raw_source_pool", "raw_build_value": int(len(ref_hash - raw_hash)), "release_reference_value": 0, "status": "若非 0，说明源快照与 release 快照不完全相同"},
    ])
    out = pd.DataFrame(rows)
    write_csv(out, output_dir / "raw_build_vs_release_comparison.csv")
    return out


def _evidence_rank(value: Any) -> int:
    text = "" if not nonempty(value) else str(value)
    return {"E3": 4, "E2": 3, "E1": 2, "E0": 1}.get(text[:2], 1)


def _system_label(value: Any) -> str:
    return SYSTEM_MAP.get("" if not nonempty(value) else str(value), "unknown/general")


def rna_length_bin(length: int) -> str:
    if length < 50:
        return "10-49"
    if length < 100:
        return "50-99"
    if length < 200:
        return "100-199"
    return "200-500"


def prediction_length_bin(length: int) -> str:
    if length < 700:
        return "short"
    if length < 1300:
        return "medium"
    return "long"


def protein_cap(size: int, rules: dict[str, Any]) -> int:
    if size <= 500:
        return rules["protein_caps"]["at_most_500"]
    if size <= 2000:
        return rules["protein_caps"]["from_501_to_2000"]
    if size <= 10000:
        return rules["protein_caps"]["from_2001_to_10000"]
    return rules["protein_caps"]["above_10000"]


def sort_for_representative(frame: pd.DataFrame) -> pd.DataFrame:
    cols = [c for c in ["evidence_rank", "has_structure", "metadata_fields", "length_deviation", "sample_id"] if c in frame]
    ascending = [False, False, False, True, True][: len(cols)]
    return frame.sort_values(cols, ascending=ascending, kind="stable")


def round_robin(frame: pd.DataFrame, strata: list[str], limit: int) -> pd.DataFrame:
    if len(frame) <= limit:
        return sort_for_representative(frame)
    buckets: list[list[Any]] = []
    actual_strata = [c for c in strata if c in frame.columns]
    if not actual_strata:
        return sort_for_representative(frame).head(limit)
    for _, group in frame.groupby(actual_strata, dropna=False, sort=True):
        buckets.append(sort_for_representative(group).index.tolist())
    selected: list[Any] = []
    cursor = 0
    while len(selected) < limit:
        progressed = False
        for bucket in buckets:
            if cursor < len(bucket):
                selected.append(bucket[cursor])
                progressed = True
                if len(selected) == limit:
                    break
        if not progressed:
            break
        cursor += 1
    return frame.loc[selected].copy()


def cap_by_protein(frame: pd.DataFrame, rules: dict[str, Any]) -> pd.DataFrame:
    chunks = []
    for _, group in frame.groupby("protein_sequence", sort=True):
        chunks.append(round_robin(group, ["rna_length_bin", "source_database"], protein_cap(len(group), rules)))
    return pd.concat(chunks, ignore_index=False) if chunks else frame.iloc[0:0].copy()


def cap_by_family(frame: pd.DataFrame, cap: int) -> pd.DataFrame:
    chunks = []
    for _, group in frame.groupby("protein_cluster30", sort=True):
        chunks.append(round_robin(group, ["rna_length_bin", "source_database"], cap))
    return pd.concat(chunks, ignore_index=False) if chunks else frame.iloc[0:0].copy()


def add_preboltz_columns(frame: pd.DataFrame, compute_expensive: bool = True, rules: dict[str, Any] | None = None) -> pd.DataFrame:
    rules = rules or load_pipeline_rules(Path(__file__).resolve().parents[1])
    out = frame.copy()
    for column, default in {
        "sample_id": "", "protein_sequence": "", "rna_sequence": "", "source_database": "",
        "protein_name": "", "system_class": "", "rna_class": "", "interaction_evidence": "",
        "experimental_structure_available": False, "benchmark_conflict": False,
        "source_filter_pass": True,
    }.items():
        if column not in out:
            out[column] = default
    out["protein_sequence"] = out["protein_sequence"].map(normalize_protein)
    out["rna_sequence"] = out["rna_sequence"].map(normalize_rna)
    out["protein_length_recalc_full"] = out["protein_sequence"].str.len().astype(int)
    out["rna_length_recalc_full"] = out["rna_sequence"].str.len().astype(int)
    if compute_expensive:
        out["recomputed_sequence_hash_full"] = [pair_hash(p, r) for p, r in zip(out["protein_sequence"], out["rna_sequence"])]
        out["hash_match_full"] = out.get("sequence_hash", "").astype(str).eq(out["recomputed_sequence_hash_full"])
    else:
        # 源加载阶段已经按归一化序列生成 sequence_hash。筛选阶段只需要它做
        # exact pair 去重；完整 hash 重算和一致性审计在最终队列中执行。
        out["recomputed_sequence_hash_full"] = out["sequence_hash"].astype(str)
        out["hash_match_full"] = True
    out["strict_sequence_pass_full"] = out["protein_sequence"].map(lambda x: bool(AA_RE.fullmatch(x))) & out["rna_sequence"].map(lambda x: bool(RNA_RE.fullmatch(x)))
    out["main_length_pass_full"] = out["protein_length_recalc_full"].between(rules["protein_length_min"], rules["protein_length_max"]) & out["rna_length_recalc_full"].between(rules["rna_length_min"], rules["rna_length_max"])
    out["evidence_level_full"] = out["interaction_evidence"].astype(str).str.extract(r"^(E[0-3])", expand=False).fillna("unknown")
    out["evidence_rank_full"] = out["interaction_evidence"].map(_evidence_rank)
    out["interaction_system_full"] = out["system_class"].map(_system_label)
    out["rna_length_bin_full"] = out["rna_length_recalc_full"].map(rna_length_bin)
    out["total_length_full"] = out["protein_length_recalc_full"] + out["rna_length_recalc_full"]
    out["prediction_length_bin_full"] = out["total_length_full"].map(prediction_length_bin)
    out["protein_invalid_chars_full"] = out["protein_sequence"].map(lambda x: "".join(sorted(set(x) - set("ACDEFGHIKLMNPQRSTVWY"))))
    out["rna_invalid_chars_full"] = out["rna_sequence"].map(lambda x: "".join(sorted(set(x) - set("ACGU"))))
    out["protein_ambiguous_count_full"] = out["protein_sequence"].map(lambda x: sum(c in AMBIGUOUS_AA for c in x))
    out["rna_ambiguous_count_full"] = out["rna_sequence"].map(lambda x: sum(c in AMBIGUOUS_RNA for c in x))
    out["rna_gc_fraction_full"] = out["rna_sequence"].map(lambda x: (x.count("G") + x.count("C")) / len(x) if x else np.nan)
    if compute_expensive:
        out["protein_entropy_full"] = out["protein_sequence"].map(entropy)
        out["rna_entropy_full"] = out["rna_sequence"].map(entropy)
    else:
        out["protein_entropy_full"] = np.nan
        out["rna_entropy_full"] = np.nan
    return out


def entropy(sequence: str) -> float:
    if not sequence:
        return float("nan")
    counts = pd.Series(list(sequence)).value_counts(normalize=True).to_numpy()
    return float(-(counts * np.log2(counts)).sum())


def load_stage(root: Path, name: str, columns: list[str] | None = None) -> pd.DataFrame:
    path = root / "data" / "final" / "pre_boltz" / name
    if not path.exists():
        raise FileNotFoundError(f"缺少阶段文件: {path}")
    return pd.read_parquet(path, columns=columns)


def load_cluster_assignments(root: Path, sequences: pd.DataFrame) -> pd.DataFrame:
    cluster_dir = root / "data" / "final" / "pre_boltz" / "clusters"
    p_path = cluster_dir / "protein_cluster_assignments.parquet"
    r_path = cluster_dir / "rna_cluster_assignments.parquet"
    if not p_path.exists() or not r_path.exists():
        raise FileNotFoundError("缺少 pre-Boltz cluster assignment 文件；请使用现有阶段快照或 --recompute-clusters。")
    p = pd.read_parquet(p_path)
    r = pd.read_parquet(r_path)
    p = p.rename(columns={"sequence": "protein_sequence", "protein90_representative": "protein_cluster90", "protein30_representative": "protein_cluster30"})
    r = r.rename(columns={"sequence": "rna_sequence", "rna90_representative": "rna_cluster90"})
    result = sequences[["protein_sequence", "rna_sequence"]].drop_duplicates().merge(p[["protein_sequence", "protein_cluster90", "protein_cluster30"]], on="protein_sequence", how="left", validate="many_to_one")
    result = result.merge(r[["rna_sequence", "rna_cluster90"]], on="rna_sequence", how="left", validate="many_to_one")
    missing = result[["protein_cluster90", "protein_cluster30", "rna_cluster90"]].isna().any(axis=1)
    if missing.any():
        raise ValueError(f"有 {int(missing.sum())} 个序列没有可复用的 pre-Boltz cluster 分配；请使用 --recompute-clusters。")
    return result


def run_mmseqs_cluster(root: Path, sequences: pd.DataFrame, mmseqs_path: Path, work_dir: Path) -> pd.DataFrame:
    """对新输入重新计算 Protein90/30 和 RNA90；参数与当前 pre_boltz.py 一致。"""
    work_dir.mkdir(parents=True, exist_ok=True)
    mmseqs = mmseqs_path.resolve()

    def one(kind: str, values: list[str], identity: float, dbtype: str) -> pd.DataFrame:
        fasta = work_dir / f"{kind}.fasta"
        mapping = {sequence: f"{kind[:1].upper()}{i:06d}" for i, sequence in enumerate(sorted(values), 1)}
        with fasta.open("w", encoding="utf-8") as handle:
            for sequence, identifier in mapping.items():
                handle.write(f">{identifier}\n{sequence}\n")
        db = work_dir / f"{kind}_db"
        clu = work_dir / f"{kind}_clu"
        tmp = work_dir / f"{kind}_tmp"
        tsv = work_dir / f"{kind}.tsv"
        commands = [
            [str(mmseqs), "createdb", str(fasta), str(db), "--dbtype", dbtype],
            [str(mmseqs), "cluster", str(db), str(clu), str(tmp), "--min-seq-id", str(identity), "-c", "0.8", "--cov-mode", "0", "--cluster-mode", "0", "--threads", "8"],
            [str(mmseqs), "createtsv", str(db), str(db), str(clu), str(tsv)],
        ]
        for command in commands:
            result = subprocess.run(command, capture_output=True, text=True)
            if result.returncode:
                detail = (result.stderr or result.stdout or "(MMseqs2 未返回错误文本)")[-4000:]
                raise RuntimeError("MMseqs2 执行失败，命令=" + " ".join(command) + "\n" + detail)
        members: dict[str, str] = {}
        for line in tsv.read_text(encoding="utf-8").splitlines():
            rep, member = line.split("\t")
            members[member] = rep
        label = f"{kind}_representative"
        return pd.DataFrame({"sequence": list(mapping), label: [members.get(mapping[s], mapping[s]) for s in mapping]})

    proteins = one("protein90", sequences["protein_sequence"].unique().tolist(), 0.90, "1")
    p30 = one("protein30", sequences["protein_sequence"].unique().tolist(), 0.30, "1")
    rnas = one("rna90", sequences["rna_sequence"].unique().tolist(), 0.90, "2")
    p = proteins.merge(p30, on="sequence", validate="one_to_one")
    p = p.rename(columns={"sequence": "protein_sequence", "protein90_representative": "protein_cluster90", "protein30_representative": "protein_cluster30"})
    r = rnas.rename(columns={"sequence": "rna_sequence", "rna90_representative": "rna_cluster90"})
    return sequences[["protein_sequence", "rna_sequence"]].drop_duplicates().merge(p, on="protein_sequence", validate="many_to_one").merge(r, on="rna_sequence", validate="many_to_one")


def rebuild_stages(root: Path, raw: pd.DataFrame, recompute_clusters: bool, mmseqs: Path | None, output_dir: Path) -> dict[str, pd.DataFrame]:
    rules = load_pipeline_rules(root)
    data = add_preboltz_columns(raw, compute_expensive=False, rules=rules)
    d1_mask = data["source_filter_pass"].fillna(True).astype(bool) & data["strict_sequence_pass_full"] & data["main_length_pass_full"] & data["sample_id"].map(nonempty) & data["source_database"].map(nonempty) & ~data["benchmark_conflict"].fillna(False).astype(bool)
    d1 = data.loc[d1_mask].copy()
    d1["hard_qc_pass"] = True
    d2 = sort_for_representative(d1.assign(evidence_rank=d1["evidence_rank_full"], has_structure=d1["experimental_structure_available"].fillna(False).astype(int), metadata_fields=0, length_deviation=0.0)).drop_duplicates("sequence_hash", keep="first").copy()
    missing_columns = [c for c in ["protein_sequence", "rna_sequence"] if c not in d2]
    if missing_columns:
        raise ValueError(f"输入表缺少列: {missing_columns}")
    if recompute_clusters:
        if mmseqs is None:
            raise ValueError("--recompute-clusters 必须同时提供 --mmseqs")
        cluster_run_dir = output_dir / "mmseqs_work" / datetime.now(timezone.utc).strftime("run_%Y%m%dT%H%M%SZ")
        mapping = run_mmseqs_cluster(root, d2[["protein_sequence", "rna_sequence"]], mmseqs, cluster_run_dir)
    else:
        mapping = load_cluster_assignments(root, d2)
    d2 = d2.drop(columns=["protein_cluster90", "protein_cluster30", "rna_cluster90"], errors="ignore").merge(mapping, on=["protein_sequence", "rna_sequence"], validate="many_to_one")
    d2["interaction_system"] = d2["system_class"].map(_system_label)
    d2["rna_length_bin"] = d2["rna_length_recalc_full"].map(rna_length_bin)
    d2["evidence_rank"] = d2["evidence_rank_full"]
    d2["has_structure"] = d2["experimental_structure_available"].fillna(False).astype(int)
    d2["metadata_fields"] = 0
    d2["local_rna_cluster_size"] = d2.groupby(["protein_sequence", "rna_cluster90"])["sample_id"].transform("size")
    medians = d2.groupby(["protein_sequence", "rna_cluster90"])["rna_length_recalc_full"].transform("median")
    d2["length_deviation"] = (d2["rna_length_recalc_full"] - medians).abs()
    d3 = sort_for_representative(d2).drop_duplicates(["protein_sequence", "rna_cluster90"], keep="first").copy()
    d3["rare_system"] = d3["interaction_system"].isin(RARE_SYSTEMS)
    rare = d3[d3["rare_system"]].copy()
    common = d3[~d3["rare_system"]].copy()
    protein_capped = cap_by_protein(common, rules)
    protein_capped_all = pd.concat([rare, protein_capped], ignore_index=False)
    family_cap = math.ceil(rules["nominal_target_rows"] * rules["protein_family_max_fraction"])
    family_capped = cap_by_family(protein_capped, family_cap)
    family_capped_all = pd.concat([rare, family_capped], ignore_index=False)
    slots = max(rules["nominal_target_rows"] - len(rare), 0)
    selected_common = round_robin(family_capped, ["interaction_system", "protein_cluster30", "rna_length_bin", "source_database"], slots)
    queue = pd.concat([rare, selected_common], ignore_index=True).drop_duplicates("sample_id").copy()
    queue["total_length"] = queue["protein_length_recalc_full"] + queue["rna_length_recalc_full"]
    queue["prediction_length_bin"] = queue["total_length"].map(prediction_length_bin)
    queue["selection_stage"] = "pre_boltz"
    queue["selection_reason"] = queue["rare_system"].map({True: "rare_system_preserved", False: "balanced_protein_family_selection"})
    return {"d1": d1, "d2": d2, "d3": d3, "protein_capped": protein_capped_all, "family_capped": family_capped_all, "queue": queue}


def stage_membership(root: Path, raw: pd.DataFrame, rebuilt: dict[str, pd.DataFrame] | None) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, int]]:
    print("[1/6] 读取筛选阶段文件...", flush=True)
    if rebuilt is None:
        stage_files = {
            "d1": "general_d1_hard_qc.parquet", "d2": "general_d2_exact_dedup.parquet",
            "d3": "general_d3_local_rna_representatives.parquet", "queue": "prediction_queue.parquet",
        }
        stages = {key: load_stage(root, value, ["sample_id"]) for key, value in stage_files.items()}
        d3_full = load_stage(root, stage_files["d3"])
        queue_full = load_stage(root, stage_files["queue"])
        # D3 已经包含 pre-Boltz 计算字段。这里只补充筛选限额需要的列，
        # 不对 143,414 条记录再次计算全量 hash/entropy。
        if "protein_length_recalc_full" not in d3_full:
            d3_full["protein_length_recalc_full"] = d3_full["protein_sequence"].astype(str).str.len()
        if "rna_length_recalc_full" not in d3_full:
            d3_full["rna_length_recalc_full"] = d3_full["rna_sequence"].astype(str).str.len()
        if "interaction_system" not in d3_full:
            d3_full["interaction_system"] = d3_full["system_class"].map(_system_label)
        if "rna_length_bin" not in d3_full:
            d3_full["rna_length_bin"] = d3_full["rna_length_recalc_full"].map(rna_length_bin)
        if "evidence_rank" not in d3_full:
            d3_full["evidence_rank"] = d3_full["interaction_evidence"].map(_evidence_rank)
        if "has_structure" not in d3_full:
            d3_full["has_structure"] = d3_full["experimental_structure_available"].fillna(False).astype(int)
        if "metadata_fields" not in d3_full:
            d3_full["metadata_fields"] = 0
        if "length_deviation" not in d3_full:
            d3_full["length_deviation"] = 0.0
        d3_full["rare_system"] = d3_full["interaction_system"].isin(RARE_SYSTEMS)
        rare = d3_full[d3_full["rare_system"]]
        common = d3_full[~d3_full["rare_system"]]
        rules = load_pipeline_rules(root)
        protein_all = pd.concat([rare, cap_by_protein(common, rules)], ignore_index=False)
        family_all = pd.concat([rare, cap_by_family(cap_by_protein(common, rules), math.ceil(rules["nominal_target_rows"] * rules["protein_family_max_fraction"]))], ignore_index=False)
        stages["protein_capped"] = protein_all[["sample_id"]]
        stages["family_capped"] = family_all[["sample_id"]]
        if "protein_length_recalc_full" not in queue_full:
            queue_full["protein_length_recalc_full"] = queue_full["protein_sequence"].astype(str).str.len()
        if "rna_length_recalc_full" not in queue_full:
            queue_full["rna_length_recalc_full"] = queue_full["rna_sequence"].astype(str).str.len()
        return _make_membership(raw, stages), queue_full, {key: len(value) for key, value in stages.items()}
    stages = {key: value for key, value in rebuilt.items()}
    membership = _make_membership(raw, {key: value[["sample_id"]] for key, value in stages.items()})
    return membership, stages["queue"], {key: len(value) for key, value in stages.items()}


def _make_membership(raw: pd.DataFrame, stages: dict[str, pd.DataFrame]) -> pd.DataFrame:
    out = raw[["sample_id"]].copy()
    for name, frame in stages.items():
        ids = set(frame["sample_id"].dropna().astype(str))
        out[f"stage_{name}"] = out["sample_id"].astype(str).isin(ids)
    ordered = ["d1", "d2", "d3", "protein_capped", "family_capped", "queue"]
    out["stage_exclusion_reason"] = ""
    for index, name in enumerate(ordered):
        failed = ~out[f"stage_{name}"]
        prior = out[[f"stage_{previous}" for previous in ordered[:index]]].all(axis=1) if index else pd.Series(True, index=out.index)
        out.loc[failed & prior & (out["stage_exclusion_reason"] == ""), "stage_exclusion_reason"] = name
    return out


def load_benchmark_sequences(root: Path) -> tuple[set[str], set[str], set[str]]:
    proteins: set[str] = set()
    rnas: set[str] = set()
    pairs: set[str] = set()
    for split in ["validation.parquet", "test.parquet", "hard_test.parquet"]:
        path = root / "data" / "reference" / "benchmark" / split
        if not path.exists():
            continue
        frame = pd.read_parquet(path, columns=["protein_sequence", "rna_sequence"])
        for protein, rna in zip(frame["protein_sequence"].map(normalize_protein), frame["rna_sequence"].map(normalize_rna)):
            proteins.add(protein)
            rnas.add(rna)
            pairs.add(pair_hash(protein, rna))
    return proteins, rnas, pairs


def enrich_final(queue: pd.DataFrame, root: Path) -> pd.DataFrame:
    # 这里必须重新计算 pair hash 做一致性核验；entropy 只是附加统计，
    # 不应让 10 万条队列的导出被不必要的逐字符计算阻塞。
    out = add_preboltz_columns(queue, compute_expensive=False, rules=load_pipeline_rules(root))
    out["recomputed_sequence_hash_full"] = [pair_hash(p, r) for p, r in zip(out["protein_sequence"], out["rna_sequence"])]
    out["hash_match_full"] = out["sequence_hash"].astype(str).eq(out["recomputed_sequence_hash_full"])
    for column in ["protein_cluster30", "protein_cluster90", "rna_cluster90", "interaction_system", "rna_length_bin", "total_length", "prediction_length_bin", "selection_stage", "selection_reason"]:
        if column not in out:
            out[column] = ""
    out["protein_pair_count_final"] = out.groupby("protein_sequence")["sample_id"].transform("size")
    out["rna_pair_count_final"] = out.groupby("rna_sequence")["sample_id"].transform("size")
    out["protein_cluster30_pair_count_final"] = out.groupby("protein_cluster30")["sample_id"].transform("size")
    out["rna_cluster90_pair_count_final"] = out.groupby("rna_cluster90")["sample_id"].transform("size")
    benchmark_proteins, benchmark_rnas, benchmark_pairs = load_benchmark_sequences(root)
    out["benchmark_exact_pair_overlap"] = out["recomputed_sequence_hash_full"].isin(benchmark_pairs)
    out["benchmark_exact_protein_overlap"] = out["protein_sequence"].isin(benchmark_proteins)
    out["benchmark_exact_rna_overlap"] = out["rna_sequence"].isin(benchmark_rnas)
    out["benchmark_conflict_recomputed"] = out["benchmark_exact_pair_overlap"] | out["benchmark_exact_protein_overlap"] | out["benchmark_exact_rna_overlap"]
    out["final_row_number"] = np.arange(1, len(out) + 1)
    return out


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, encoding="utf-8-sig")


def write_analysis(
    root: Path,
    output_dir: Path,
    raw: pd.DataFrame,
    final: pd.DataFrame,
    membership: pd.DataFrame,
    stage_counts: dict[str, int],
    download_results: list[dict[str, Any]] | None,
    input_path: str,
    build_origin: str,
    reference_comparison: pd.DataFrame | None = None,
) -> dict[str, Any]:
    rules = load_pipeline_rules(root)
    output_dir.mkdir(parents=True, exist_ok=True)
    waterfall = pd.DataFrame([
        {"stage": "input_general", "rows": len(raw), "retention_vs_input": 1.0},
        {"stage": "D1_hard_qc", "rows": stage_counts.get("d1", 0), "retention_vs_input": stage_counts.get("d1", 0) / max(len(raw), 1)},
        {"stage": "D2_exact_dedup", "rows": stage_counts.get("d2", 0), "retention_vs_input": stage_counts.get("d2", 0) / max(len(raw), 1)},
        {"stage": "D3_local_RNA90_representatives", "rows": stage_counts.get("d3", 0), "retention_vs_input": stage_counts.get("d3", 0) / max(len(raw), 1)},
        {"stage": "protein_caps", "rows": stage_counts.get("protein_capped", 0), "retention_vs_input": stage_counts.get("protein_capped", 0) / max(len(raw), 1)},
        {"stage": "family_caps", "rows": stage_counts.get("family_capped", 0), "retention_vs_input": stage_counts.get("family_capped", 0) / max(len(raw), 1)},
        {"stage": "prediction_queue", "rows": len(final), "retention_vs_input": len(final) / max(len(raw), 1)},
    ])
    write_csv(waterfall, output_dir / "stage_waterfall.csv")

    source = final.assign(source_database=final["source_database"].fillna("<missing>")).groupby("source_database").size().rename("pairs").reset_index()
    source["fraction"] = source["pairs"] / max(len(final), 1)
    write_csv(source.sort_values("pairs", ascending=False), output_dir / "source_distribution.csv")

    protein_degree = final.groupby("protein_sequence", as_index=False).agg(
        pair_count=("sample_id", "size"),
        protein_name=("protein_name", lambda values: next((str(v) for v in values if nonempty(v)), "")),
        protein_cluster30=("protein_cluster30", "first"),
        protein_cluster90=("protein_cluster90", "first"),
        protein_length=("protein_length_recalc_full", "first"),
    ).sort_values(["pair_count", "protein_name"], ascending=[False, True])
    protein_degree.insert(0, "rank", np.arange(1, len(protein_degree) + 1))
    write_csv(protein_degree, output_dir / "protein_degree.csv")

    rna_degree = final.groupby("rna_sequence", as_index=False).agg(
        pair_count=("sample_id", "size"),
        rna_class=("rna_class", "first"),
        rna_cluster90=("rna_cluster90", "first"),
        rna_length=("rna_length_recalc_full", "first"),
    ).sort_values("pair_count", ascending=False)
    rna_degree.insert(0, "rank", np.arange(1, len(rna_degree) + 1))
    write_csv(rna_degree, output_dir / "rna_degree.csv")

    raw_columns = list(raw.columns)
    completeness_rows = []
    for column in raw_columns:
        series = raw[column]
        null_count = int(series.isna().sum())
        blank_count = int((~series.isna() & series.astype(str).str.strip().eq("")).sum())
        completeness_rows.append({"field": column, "rows": len(raw), "null_count": null_count, "blank_string_count": blank_count, "nonempty_count": len(raw) - null_count - blank_count, "nonempty_fraction": (len(raw) - null_count - blank_count) / max(len(raw), 1)})
    write_csv(pd.DataFrame(completeness_rows), output_dir / "field_completeness.csv")

    quality = pd.DataFrame([
        {"metric": "final_rows", "value": len(final), "definition": "最终 pre-Boltz 队列行数"},
        {"metric": "unique_proteins", "value": final["protein_sequence"].nunique(), "definition": "最终队列中的不同蛋白序列数"},
        {"metric": "unique_rnas", "value": final["rna_sequence"].nunique(), "definition": "最终队列中的不同 RNA 序列数"},
        {"metric": "unique_pairs", "value": final["recomputed_sequence_hash_full"].nunique(), "definition": "蛋白序列与 RNA 序列联合 hash 的不同值数量"},
        {"metric": "hash_mismatch_rows", "value": int((~final["hash_match_full"]).sum()), "definition": "sequence_hash 与当前序列重算 hash 不一致的行数"},
        {"metric": "duplicate_sample_id_rows", "value": int(final["sample_id"].duplicated().sum()), "definition": "最终队列重复 sample_id 行数"},
        {"metric": "benchmark_exact_pair_overlap_rows", "value": int(final["benchmark_exact_pair_overlap"].sum()), "definition": "与 benchmark 同时完全相同 protein+RNA pair 的行数"},
        {"metric": "benchmark_exact_protein_overlap_rows", "value": int(final["benchmark_exact_protein_overlap"].sum()), "definition": "与 benchmark 完全相同蛋白序列的行数"},
        {"metric": "benchmark_exact_rna_overlap_rows", "value": int(final["benchmark_exact_rna_overlap"].sum()), "definition": "与 benchmark 完全相同 RNA 序列的行数"},
        {"metric": "strict_sequence_pass_rows", "value": int(final["strict_sequence_pass_full"].sum()), "definition": "蛋白只含 20 种标准氨基酸且 RNA 只含 A/C/G/U"},
        {"metric": "main_length_pass_rows", "value": int(final["main_length_pass_full"].sum()), "definition": f"蛋白 {rules['protein_length_min']}–{rules['protein_length_max']} aa 且 RNA {rules['rna_length_min']}–{rules['rna_length_max']} nt"},
        {"metric": "rfam_annotated_rows", "value": int(final.get("rfam_family", pd.Series(index=final.index, dtype=object)).map(nonempty).sum()), "definition": "有 Rfam 注释的最终行数"},
    ])
    write_csv(quality, output_dir / "quality_summary.csv")

    stage_audit = raw.merge(membership, on="sample_id", how="left", validate="one_to_one")
    write_csv(stage_audit, output_dir / "general_stage_audit.csv")

    data_dictionary = pd.DataFrame({
        "field": list(final.columns),
        "meaning": [
            "最终队列保留的原始字段或 pre-Boltz 字段；详细口径见脚本源代码。" if not c.endswith("_full") else "本脚本由当前序列重新计算的审计字段。"
            for c in final.columns
        ],
    })
    write_csv(data_dictionary, output_dir / "final_csv_data_dictionary.csv")

    pdb_ids = final.get("pdb_id", pd.Series(index=final.index, dtype=object)).fillna("").astype(str).str.strip()
    experimental_structures = int(pdb_ids.replace("", pd.NA).nunique())
    predicted_complexes = int(final.loc[pdb_ids.eq(""), "recomputed_sequence_hash_full"].nunique())
    summary = {
        "script_version": SCRIPT_VERSION,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "input_path": input_path,
        "build_origin": build_origin,
        "output_dir": str(output_dir.relative_to(root)),
        "raw_source_candidate_rows": int(len(raw)),
        "final_csv_rows": int(len(final)),
        "final_csv_columns": int(len(final.columns)),
        "experimental_pdb_structures": experimental_structures,
        "predicted_complexes": predicted_complexes,
        "dataset_scope": "single_PRI-General_table_with_experimental_PDB_and_predicted_complex_units",
        "stage_counts": {key: int(value) for key, value in stage_counts.items()},
        "unique_final_proteins": int(final["protein_sequence"].nunique()),
        "unique_final_rnas": int(final["rna_sequence"].nunique()),
        "unique_final_pairs": int(final["recomputed_sequence_hash_full"].nunique()),
        "source_distribution": {str(k): int(v) for k, v in final["source_database"].fillna("<missing>").value_counts().items()},
        "evidence_distribution": {str(k): int(v) for k, v in final["evidence_level_full"].value_counts().items()},
        "rna_class_distribution": {str(k): int(v) for k, v in final["rna_class"].fillna("<missing>").value_counts().items()},
        "benchmark_exact_overlap": {
            "pair": int(final["benchmark_exact_pair_overlap"].sum()),
            "protein": int(final["benchmark_exact_protein_overlap"].sum()),
            "rna": int(final["benchmark_exact_rna_overlap"].sum()),
        },
        "download_mode": download_results is not None,
        "download_results": download_results or [],
        "rules": {
            "protein_length": f"{rules['protein_length_min']}-{rules['protein_length_max']} aa",
            "rna_length": f"{rules['rna_length_min']}-{rules['rna_length_max']} nt",
            "protein_alphabet": "20 canonical amino acids only",
            "rna_alphabet": "A/C/G/U only",
            "protein_caps": rules["protein_caps"],
            "protein30_family_max_fraction": rules["protein_family_max_fraction"],
            "nominal_target_rows": rules["nominal_target_rows"],
            "family_cap_rows": math.ceil(rules["nominal_target_rows"] * rules["protein_family_max_fraction"]),
            "boltz_status": "not run",
        },
    }
    if reference_comparison is not None:
        summary["release_reference_comparison"] = {
            str(row["metric"]): {
                "raw_build_value": row.get("raw_build_value"),
                "release_reference_value": row.get("release_reference_value"),
                "status": row.get("status"),
            }
            for _, row in reference_comparison.iterrows()
        }
    (output_dir / "pipeline_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    report = [
        "# PRI-General 单文件流程输出",
        "",
        "当前项目只发布一个 PRI-General 数据集，即本次脚本生成的最终表。表内包含两类结构单位：已有 PDB 实验结构，以及没有实验结构、后续需要预测的复合物。",
        "",
        f"已有 PDB 实验结构：**{experimental_structures:,} 个**；待预测复合物：**{predicted_complexes:,} 个**。最终表中的技术性记录数为 {len(final):,}，不作为结构数使用。",
        "",
        f"本次构建来源：`{build_origin}`。未运行 Boltz，也未修改 data/final 中的最终快照。Benchmark 只作为训练排除和泄漏检查的参考集合，不属于 PRI-General 发布数据。",
        "",
        "## 输出文件",
        "",
        "- `PRI-General_pre_boltz_detailed.csv`：最终队列，每行一条 pair，保留原始字段并附加长度、GC、hash、质量、cluster、degree 和 benchmark 审计字段。",
        "- `general_stage_audit.csv`：源级候选记录的逐阶段进入/未进入标记。",
        "- `source_ingestion_summary.csv`：四类来源实际读取行数及来源级过滤结果。",
        "- `raw_build_vs_release_comparison.csv`：raw-to-final 结果与 release General 的独立对照，不作为输入。",
        "- `stage_waterfall.csv`：从输入到最终队列的数量瀑布。",
        "- `protein_degree.csv`、`rna_degree.csv`：最终队列中每个唯一蛋白/RNA 的配对次数。",
        "- `field_completeness.csv`、`quality_summary.csv`：字段完整性和质量统计。",
        "",
        "## 运行方式",
        "",
        "`python scripts/pri_general_full_pipeline.py`",
        "",
        "默认从 `data/source_candidates/` 开始重建；若要下载可直接下载的来源快照，在命令后加入 `--download`。如需重新聚类，请显式提供目标环境中的 MMseqs2。",
    ]
    (output_dir / "report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="PRI-General 单文件构建、审计和完整 CSV 导出")
    parser.add_argument("--root", type=Path, default=None, help="项目根目录，默认取脚本上两级目录")
    parser.add_argument("--input", type=Path, default=None, help="复核模式下使用的标准化 General parquet")
    parser.add_argument("--raw-candidates-dir", type=Path, default=None, help="源级候选 parquet 目录；默认 data/source_candidates")
    parser.add_argument("--from-release-snapshot", action="store_true", help="仅复核模式：把 release General 作为输入；默认不启用")
    parser.add_argument("--output-dir", type=Path, default=None, help="分析输出目录")
    parser.add_argument("--final-csv", type=Path, default=None, help="最终详细 CSV 的路径")
    parser.add_argument("--download", action="store_true", help="按脚本内 SOURCE_DOWNLOAD_PLAN 下载可直接下载的来源")
    parser.add_argument("--download-dir", type=Path, default=None, help="下载目录")
    parser.add_argument("--force-download", action="store_true", help="覆盖已有下载文件")
    parser.add_argument("--rebuild", action="store_true", help="复核模式下重新执行 D1-D3、蛋白限额和家族限额；源数据模式自动执行")
    parser.add_argument("--recompute-clusters", action="store_true", help="使用 MMseqs2 重新计算聚类")
    parser.add_argument("--mmseqs", type=Path, default=None, help="MMseqs2 可执行文件")
    args = parser.parse_args()

    root = project_root(args.root)
    output_dir = (args.output_dir or root / "reports" / "final" / "rebuild").resolve()
    final_csv = (args.final_csv or root / "data" / "final" / "rebuild" / "PRI-General_pre_boltz_detailed.csv").resolve()

    download_results = download_sources(root, args.download_dir, args.force_download) if args.download else None
    raw_source_summary = None
    if args.from_release_snapshot:
        input_path = (args.input or root / "data" / "reference" / "general_source_snapshot.parquet").resolve()
        if not input_path.exists():
            raise FileNotFoundError(f"输入文件不存在: {input_path}")
        print(f"[0/6] 复核模式读取 release 快照: {input_path}", flush=True)
        raw = pd.read_parquet(input_path)
        build_origin = "release_snapshot_reference_only"
    else:
        candidates_dir = (args.raw_candidates_dir or root / "data" / "source_candidates").resolve()
        print(f"[0/6] 从源级候选开始构建: {candidates_dir}", flush=True)
        raw, raw_source_summary = load_raw_source_candidates(root, candidates_dir, output_dir)
        input_path = candidates_dir
        build_origin = "raw_source_candidates"
        benchmark_proteins, benchmark_rnas, benchmark_pairs = load_benchmark_sequences(root)
        raw["benchmark_conflict"] = [
            pair_hash(p, r) in benchmark_pairs or p in benchmark_proteins or r in benchmark_rnas
            for p, r in zip(raw["protein_sequence"], raw["rna_sequence"])
        ]
        print(f"[raw] benchmark 冲突候选: {int(raw['benchmark_conflict'].sum()):,} 条；将在 D1 前排除", flush=True)
    if "sample_id" not in raw or raw["sample_id"].duplicated().any():
        raise ValueError("输入 General 表必须有唯一 sample_id。")
    required = {"protein_sequence", "rna_sequence", "source_database"}
    missing = sorted(required - set(raw.columns))
    if missing:
        raise ValueError("输入 General 表缺少列: " + ", ".join(missing))

    rebuilt = None
    if not args.from_release_snapshot or args.rebuild:
        print("[1/6] 重新执行 D1-D3 和限额筛选...", flush=True)
        mmseqs = args.mmseqs
        recompute = args.recompute_clusters
        if not args.from_release_snapshot and mmseqs is None:
            # 独立工作区不捆绑 MMseqs2 二进制；没有显式提供时复用当前
            # 工作区中的最终 cluster 分配。
            recompute = False
            print("[cluster] 未提供 MMseqs2，复用工作区内已有 cluster 分配。", flush=True)
        rebuilt = rebuild_stages(root, raw, recompute, mmseqs, output_dir)
    membership, queue, stage_counts = stage_membership(root, raw, rebuilt)
    print(f"[2/6] 阶段完成: input={len(raw):,}, D3={stage_counts.get('d3', 0):,}, queue={len(queue):,}", flush=True)
    print("[3/6] 计算最终队列审计字段...", flush=True)
    final = enrich_final(queue, root)
    print(f"[4/6] 写出最终 CSV: {final_csv}（{len(final):,} 条）", flush=True)
    write_csv(final, final_csv)
    print("[5/6] 生成统计表和阶段审计表...", flush=True)
    reference_comparison = compare_with_release_general(root, raw, final, output_dir) if not args.from_release_snapshot else None
    summary = write_analysis(root, output_dir, raw, final, membership, stage_counts, download_results, str(input_path), build_origin, reference_comparison)
    print("[6/6] 完成。", flush=True)
    print(json.dumps({"final_csv": str(final_csv), "analysis_dir": str(output_dir), **summary}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
