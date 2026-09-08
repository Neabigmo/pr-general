"""Run the four-source proposal as a sequence-resolved audit pass.

This is deliberately separate from the released v1 rebuild.  The current
workspace has sequence-ready PDB/NPInter tables, while the local RNAInter,
IntAct, and ENCODE/RBNS assets are source material without the two exact
sequences required by the proposal.  Those assets are inventoried and held;
they are never silently converted into candidate pairs.
"""

from __future__ import annotations

import argparse
import json
import re
import tarfile
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd


SOURCE_FILES = {
    "npinter_main": "npinter_raw.parquet",
    "npinter_bindingsite": "npinter_bindingsite_raw.parquet",
    "npinter_mirna": "mirna_protein_raw.parquet",
    "rcsb_pdb": "pdb_raw.parquet",
}

AA_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")
RNA_RE = re.compile(r"^[ACGU]+$")
PUBMED_RE = re.compile(r"(?i)(?:pubmed|pmid|reference|ref)\s*[=:]\s*(\d{6,9})")


def clean(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    return "" if text.lower() in {"nan", "none", "null"} else text


def source_category(source_database: str) -> str:
    return "PDB" if source_database == "RCSB_PDB" else "NPInter"


def evidence_rank(row: pd.Series) -> int:
    if bool(row["experimental_structure"]):
        return 4
    text = f"{row['interaction_evidence']} {row['notes']}".lower()
    if text.startswith("e3") or "e3:" in text:
        return 3
    if text.startswith("e2") or "e2:" in text:
        return 2
    return 1


def assay_type(row: pd.Series) -> str:
    if bool(row["experimental_structure"]):
        return "experimental_structure"
    text = f"{row['interaction_evidence']} {row['notes']}".lower()
    if "binding-site" in text or "bindingsite" in text or "clip" in text:
        return "CLIP/site-window"
    if "mirna" in text:
        return "miRNA-protein pairing"
    if any(token in text for token in ("emsa", "immunoprecipitation", "biochemical")):
        return "biochemical/curated"
    return "source-reported interaction"


def publication_id(row: pd.Series) -> str:
    match = PUBMED_RE.search(f"{row['interaction_evidence']} {row['notes']}")
    return match.group(1) if match else ""


def load_candidates(root: Path) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for source_id, filename in SOURCE_FILES.items():
        path = root / "data" / "source_candidates" / filename
        if not path.exists():
            raise FileNotFoundError(path)
        frame = pd.read_parquet(path)
        frame["input_source_id"] = source_id
        frame["input_file"] = filename
        frames.append(frame)
    raw = pd.concat(frames, ignore_index=True)
    raw["protein_sequence"] = raw["protein_sequence"].map(clean).str.upper()
    raw["rna_sequence"] = raw["rna_sequence"].map(clean).str.upper().str.replace("T", "U", regex=False)
    raw["protein_id"] = raw["protein_accession"].map(clean)
    raw["rna_id"] = raw["rna_accession"].map(clean)
    raw["species"] = raw["organism"].map(clean)
    raw["source_database"] = raw["source_database"].map(clean)
    raw["source_record_id"] = raw["source_record_id"].map(clean)
    raw["pdb_id"] = raw["pdb_id"].map(clean)
    raw["experimental_structure"] = raw["experimental_structure_available"].fillna(False).astype(bool)
    raw["proposal_source"] = raw["source_database"].map(source_category)
    raw["exact_protein_mapping"] = raw["protein_sequence"].ne("") & raw["protein_id"].ne("")
    raw["exact_rna_mapping"] = raw["rna_sequence"].ne("") & raw["rna_id"].ne("")
    raw["protein_length"] = raw["protein_sequence"].str.len()
    raw["rna_length"] = raw["rna_sequence"].str.len()
    raw["valid_alphabet"] = raw["protein_sequence"].map(lambda value: bool(AA_RE.fullmatch(value))) & raw["rna_sequence"].map(lambda value: bool(RNA_RE.fullmatch(value)))
    raw["length_pass"] = raw["protein_length"].between(40, 2000) & raw["rna_length"].between(12, 500)
    raw["pdb_contact_pass"] = (~raw["proposal_source"].eq("PDB")) | raw["interaction_evidence"].fillna("").str.contains("contact atom pairs", case=False, regex=False)
    raw["pdb_record_pass"] = (~raw["proposal_source"].eq("PDB")) | (raw["experimental_structure"] & raw["pdb_id"].ne(""))
    raw["candidate_eligible"] = raw["exact_protein_mapping"] & raw["exact_rna_mapping"] & raw["valid_alphabet"] & raw["length_pass"] & raw["pdb_contact_pass"] & raw["pdb_record_pass"]
    raw["hold_reason"] = ""
    raw.loc[~raw["exact_protein_mapping"], "hold_reason"] = "missing_exact_protein_mapping"
    raw.loc[raw["hold_reason"].eq("") & ~raw["exact_rna_mapping"], "hold_reason"] = "missing_exact_rna_mapping"
    raw.loc[raw["hold_reason"].eq("") & ~raw["valid_alphabet"], "hold_reason"] = "invalid_sequence_alphabet"
    raw.loc[raw["hold_reason"].eq("") & ~raw["length_pass"], "hold_reason"] = "proposal_length_filter"
    raw.loc[raw["hold_reason"].eq("") & ~raw["pdb_contact_pass"], "hold_reason"] = "pdb_contact_evidence_missing"
    raw.loc[raw["hold_reason"].eq("") & ~raw["pdb_record_pass"], "hold_reason"] = "pdb_record_not_experimental"
    raw["evidence_type"] = np.where(raw["experimental_structure"], "experimental_structure", raw["interaction_evidence"].fillna("").str.extract(r"^(E\d+)", expand=False).fillna("source_reported"))
    raw["assay_type"] = raw.apply(assay_type, axis=1)
    raw["publication_id"] = raw.apply(publication_id, axis=1)
    raw["pair_key"] = raw["protein_sequence"] + "\x1f" + raw["rna_sequence"]
    return raw


def make_pair_pool(records: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    eligible = records.loc[records["candidate_eligible"]].copy()
    eligible["evidence_rank"] = eligible.apply(evidence_rank, axis=1)
    eligible = eligible.sort_values(
        ["pair_key", "evidence_rank", "experimental_structure", "source_record_id"],
        ascending=[True, False, False, True],
    )
    evidence = eligible.copy()
    evidence["source_category"] = evidence["proposal_source"]
    pair_counts = evidence.groupby("pair_key", sort=False).agg(
        evidence_record_count=("source_record_id", "size"),
        source_count=("proposal_source", "nunique"),
        source_categories=("proposal_source", lambda values: ";".join(sorted(set(values)))),
        pdb_ids=("pdb_id", lambda values: ";".join(sorted({value for value in values if value}))),
    )
    pairs = evidence.drop_duplicates("pair_key", keep="first").merge(pair_counts, left_on="pair_key", right_index=True, how="left", validate="one_to_one")
    pairs["source_category"] = pairs["proposal_source"]
    pairs["source_category"] = pairs["source_category"].map(clean)
    return pairs, evidence


def source_inventory(root: Path, external_root: Path, records: pd.DataFrame) -> list[dict[str, object]]:
    current = root / "data" / "source_candidates"
    eclip_root = root / "data" / "downloads" / "encode_eclip"
    current_counts = records["proposal_source"].value_counts().to_dict()
    eclip_status = "site_asset_ready" if (eclip_root / "manifest.json").exists() else "missing"
    inventory: list[dict[str, object]] = []
    for label, path, kind, status, note in [
        ("PDB", current / SOURCE_FILES["rcsb_pdb"], "sequence/structure candidate table", "candidate_ready", "current source table"),
        ("NPInter", current / SOURCE_FILES["npinter_main"], "sequence candidate table", "candidate_ready", "main table; unresolved rows are held by exact mapping filter"),
        ("RNAInter", external_root / "rnainter", "raw interaction archives", "hold", "symbol/ID relations; no exact double-sequence canonical table"),
        ("IntAct", external_root / "intact" / "intact.zip", "MITAB interaction archive", "hold", "RNA sequence resolution is not present in the local MITAB asset"),
        ("ENCODE/RBNS", external_root / "encode_rbns", "RBNS/RNAcompete metadata and matrices", "assay_only", "local asset is not an eCLIP genomic site-window table"),
        ("ENCODE/eCLIP", eclip_root, "reproducible peak BED files", eclip_status, "requires genome+transcript+strand window extraction and unique RBP protein mapping"),
    ]:
        exists = path.exists()
        row: dict[str, object] = {
            "source_category": label,
            "local_path": str(path.resolve()),
            "exists": exists,
            "asset_kind": kind,
            "candidate_pool_status": status,
            "exact_mapping_status": "sequence_resolved" if status == "candidate_ready" else "not_sequence_resolved",
            "observed_assets": 1 if exists and path.is_file() else 0,
            "observed_rows": None,
            "note": note,
        }
        if label in current_counts:
            row["observed_rows"] = int(current_counts[label])
        if label == "RNAInter" and path.is_dir():
            files = sorted(path.glob("*.tar.gz"))
            row["observed_assets"] = len(files)
            row["observed_rows"] = sum(count_tar_rows(file) for file in files)
        elif label == "IntAct" and exists:
            with zipfile.ZipFile(path) as archive:
                row["observed_assets"] = len(archive.infolist())
                row["observed_rows"] = legacy_hold_rows(external_root)
        elif label == "ENCODE/RBNS" and path.is_dir():
            metadata = path / "file_metadata"
            processed = path / "processed"
            row["observed_assets"] = len(list(metadata.glob("*.json"))) + len(list(processed.glob("*")))
            row["observed_rows"] = len(list(processed.glob("*")))
        elif label == "ENCODE/eCLIP" and path.is_dir():
            manifest_path = path / "manifest.json"
            if manifest_path.exists():
                manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
                row["observed_assets"] = len(manifest.get("files", []))
                row["observed_rows"] = len(manifest.get("files", []))
                row["exact_mapping_status"] = "genomic_sites_not_sequence_resolved"
        inventory.append(row)
    return inventory


def count_tar_rows(path: Path) -> int:
    with tarfile.open(path, mode="r:gz") as archive:
        member = archive.next()
        if member is None:
            return 0
        handle = archive.extractfile(member)
        if handle is None:
            return 0
        newlines = 0
        while chunk := handle.read(1024 * 1024):
            newlines += chunk.count(b"\n")
    return max(newlines - 1, 0)


def legacy_hold_rows(external_root: Path) -> int | None:
    hold_candidates = [
        external_root.parent.parent.parent / "1DATA" / "data" / "hold" / "intact.parquet",
        Path(r"H:\2026try\8.15PR\1DATA\data\hold\intact.parquet"),
    ]
    for path in hold_candidates:
        if path.exists():
            return int(len(pd.read_parquet(path, columns=[])))
    return None


def stage_summary(records: pd.DataFrame, pairs: pd.DataFrame) -> pd.DataFrame:
    stages = [
        ("raw_records", len(records), "all four current sequence tables"),
        ("exact_protein_mapping", int(records["exact_protein_mapping"].sum()), "protein sequence and identifier both present"),
        ("exact_rna_mapping", int(records["exact_rna_mapping"].sum()), "RNA sequence and identifier both present"),
        ("exact_double_mapping", int((records["exact_protein_mapping"] & records["exact_rna_mapping"]).sum()), "both exact mapping booleans true"),
        ("alphabet_and_length_qc", int((records["exact_protein_mapping"] & records["exact_rna_mapping"] & records["valid_alphabet"] & records["length_pass"]).sum()), "proposal sequence alphabet and length thresholds"),
        ("pdb_contact_and_structure_qc", int(records["candidate_eligible"].sum()), "plus PDB contact/experimental-record checks"),
        ("unique_exact_sequence_pairs", int(len(pairs)), "one representative; evidence records retained separately"),
    ]
    return pd.DataFrame(stages, columns=["stage", "rows", "definition"])


def source_funnel(records: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for source_id, group in records.groupby("input_source_id", sort=True):
        double = group["exact_protein_mapping"] & group["exact_rna_mapping"]
        eligible = group["candidate_eligible"]
        protein_alphabet = group["protein_sequence"].map(lambda value: bool(re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", value)))
        rna_alphabet = group["rna_sequence"].map(lambda value: bool(re.fullmatch(r"[ACGU]+", value)))
        protein_length = group["protein_length"].between(40, 2000)
        rna_length = group["rna_length"].between(12, 500)
        rows.append({
            "input_source_id": source_id,
            "raw_rows": len(group),
            "protein_exact_rows": int(group["exact_protein_mapping"].sum()),
            "rna_exact_rows": int(group["exact_rna_mapping"].sum()),
            "double_exact_rows": int(double.sum()),
            "protein_alphabet_pass_rows": int((double & protein_alphabet).sum()),
            "rna_alphabet_pass_rows": int((double & rna_alphabet).sum()),
            "protein_length_pass_rows": int((double & protein_length).sum()),
            "rna_length_pass_rows": int((double & rna_length).sum()),
            "alphabet_length_pass_rows": int((double & group["valid_alphabet"] & group["length_pass"]).sum()),
            "eligible_evidence_rows": int(eligible.sum()),
            "unique_eligible_pairs": int(group.loc[eligible, "pair_key"].nunique()),
        })
    return pd.DataFrame(rows)


def npinter_resolution_profile(records: pd.DataFrame) -> pd.DataFrame:
    main = records.loc[records["input_source_id"].eq("npinter_main")].copy()
    main["protein_resolution_source"] = main["notes"].fillna("").str.extract(r"protein=([^;]+)", expand=False).fillna("<unlabeled>").map(clean)
    main["rna_resolution_source"] = main["notes"].fillna("").str.extract(r"RNA=([^;]+)", expand=False).fillna("<unlabeled>").map(clean)
    main["protein_id_type"] = np.select(
        [main["protein_id"].str.match(r"^(P|Q|O|A0A)", case=False), main["protein_id"].str.match(r"^(NP|XP|YP|WP)_", case=False), main["protein_id"].str.match(r"^ENS", case=False)],
        ["uniprot_like", "refseq_protein", "ensembl_like"],
        default="other_or_missing",
    )
    main["rna_id_type"] = np.select(
        [main["rna_id"].str.match(r"^NON", case=False), main["rna_id"].str.match(r"^ENS", case=False), main["rna_id"].str.match(r"^(MI|[A-Z]{3}-mir)", case=False)],
        ["noncode_like", "ensembl_like", "mirbase_like"],
        default="other_or_missing",
    )
    main["protein_mapping_state"] = np.select(
        [main["exact_protein_mapping"], main["protein_sequence"].ne(""), main["protein_id"].ne("")],
        ["exact", "sequence_without_id", "id_without_sequence"],
        default="missing",
    )
    main["rna_mapping_state"] = np.select(
        [main["exact_rna_mapping"], main["rna_sequence"].ne(""), main["rna_id"].ne("")],
        ["exact", "sequence_without_id", "id_without_sequence"],
        default="missing",
    )
    return (
        main.groupby(["protein_resolution_source", "rna_resolution_source", "protein_id_type", "rna_id_type", "protein_mapping_state", "rna_mapping_state"], dropna=False)
        .size()
        .rename("rows")
        .reset_index()
        .sort_values("rows", ascending=False)
    )


def cap_audit(pairs: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    total = len(pairs)
    source_counts = pairs["source_category"].value_counts()
    for source, count in source_counts.items():
        rows.append({"metric": f"source:{source}", "observed": float(count / max(total, 1)), "observed_rows": int(count), "limit": 0.30, "status": "OBSERVE", "note": "one source upper bound; temporarily not a hard gate"})
    max_per_protein = int(pairs.groupby("protein_sequence").size().max()) if len(pairs) else 0
    over_proteins = int((pairs.groupby("protein_sequence").size() > 250).sum()) if len(pairs) else 0
    rows.append({"metric": "exact_protein_pair_degree", "observed": max_per_protein, "observed_rows": over_proteins, "limit": 250, "status": "PASS" if max_per_protein <= 250 else "FAIL", "note": "maximum exact-sequence pairs per protein"})
    clip = pairs["assay_type"].eq("CLIP/site-window")
    rows.append({"metric": "clip_site_window_fraction", "observed": float(clip.mean()) if len(pairs) else 0.0, "observed_rows": int(clip.sum()), "limit": 0.35, "status": "PASS" if (float(clip.mean()) if len(pairs) else 0.0) <= 0.35 else "FAIL", "note": "CLIP/site-window upper bound"})
    class_fraction = pairs["rna_class"].fillna("<missing>").map(clean).replace("", "<missing>").value_counts(normalize=True)
    top_class = class_fraction.index[0] if len(class_fraction) else ""
    rows.append({"metric": "single_rna_class_fraction", "observed": float(class_fraction.iloc[0]) if len(class_fraction) else 0.0, "observed_rows": int(pairs["rna_class"].fillna("").eq(top_class).sum()) if top_class else 0, "limit": 0.30, "status": "PASS" if (float(class_fraction.iloc[0]) if len(class_fraction) else 0.0) <= 0.30 else "FAIL", "note": top_class})
    species_fraction = pairs["species"].fillna("<missing>").map(clean).replace("", "<missing>").value_counts(normalize=True)
    top_species = species_fraction.index[0] if len(species_fraction) else ""
    rows.append({"metric": "single_species_fraction", "observed": float(species_fraction.iloc[0]) if len(species_fraction) else 0.0, "observed_rows": int(pairs["species"].fillna("").eq(top_species).sum()) if top_species else 0, "limit": 0.55, "status": "PASS" if (float(species_fraction.iloc[0]) if len(species_fraction) else 0.0) <= 0.55 else "FAIL", "note": top_species})
    human_mouse = pairs["species"].fillna("").str.contains("Homo sapiens|Mus musculus", case=False, regex=True)
    rows.append({"metric": "human_mouse_fraction", "observed": float(human_mouse.mean()) if len(pairs) else 0.0, "observed_rows": int(human_mouse.sum()), "limit": 0.75, "status": "PASS" if (float(human_mouse.mean()) if len(pairs) else 0.0) <= 0.75 else "FAIL", "note": "combined human and mouse"})
    if len(source_counts) >= 2:
        upper_bound = int(min(count / 0.30 for count in source_counts.tolist()))
    else:
        upper_bound = 0
    rows.append({"metric": "source_cap_pool_upper_bound", "observed": upper_bound, "observed_rows": upper_bound, "limit": 100000, "status": "OBSERVE", "note": "diagnostic only until RNAInter/IntAct/eCLIP yield is known"})
    rows.append({"metric": "strict_family_split", "observed": None, "observed_rows": None, "limit": "required", "status": "BLOCKED", "note": "current candidate tables do not contain Protein40/RNA80/Rfam assignments"})
    return pd.DataFrame(rows)


def markdown_table(frame: pd.DataFrame) -> str:
    columns = [str(column) for column in frame.columns]
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join("---" for _ in columns) + " |"]
    for values in frame.itertuples(index=False, name=None):
        cells = ["" if value is None or (isinstance(value, float) and np.isnan(value)) else str(value).replace("|", "\\|") for value in values]
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def write_report(output: Path, records: pd.DataFrame, pairs: pd.DataFrame, stages: pd.DataFrame, funnel: pd.DataFrame, resolution_profile: pd.DataFrame, inventory: pd.DataFrame, caps: pd.DataFrame) -> dict[str, object]:
    output.mkdir(parents=True, exist_ok=True)
    pair_columns = [
        "protein_sequence", "rna_sequence", "protein_id", "rna_id", "species", "rna_class", "evidence_type", "assay_type",
        "proposal_source", "source_database", "source_record_id", "publication_id", "exact_protein_mapping", "exact_rna_mapping",
        "experimental_structure", "pdb_id", "protein_length", "rna_length", "evidence_record_count", "source_count", "source_categories", "pdb_ids",
    ]
    pairs[pair_columns].to_parquet(output / "proposal_exact_pair_pool.parquet", index=False)
    evidence_columns = [
        "protein_sequence", "rna_sequence", "protein_id", "rna_id", "species", "rna_class", "evidence_type", "assay_type",
        "proposal_source", "source_database", "source_record_id", "publication_id", "exact_protein_mapping", "exact_rna_mapping",
        "experimental_structure", "pdb_id", "protein_length", "rna_length", "input_source_id", "input_file",
    ]
    records.loc[records["candidate_eligible"], evidence_columns].to_parquet(output / "proposal_evidence_records.parquet", index=False)
    stages.to_csv(output / "proposal_stage_summary.csv", index=False, encoding="utf-8-sig")
    funnel.to_csv(output / "proposal_source_funnel.csv", index=False, encoding="utf-8-sig")
    resolution_profile.to_csv(output / "npinter_main_resolution_profile.csv", index=False, encoding="utf-8-sig")
    inventory.to_csv(output / "proposal_source_readiness.csv", index=False, encoding="utf-8-sig")
    caps.to_csv(output / "proposal_cap_audit.csv", index=False, encoding="utf-8-sig")
    summary = {
        "run_type": "four_source_proposal_audit",
        "current_sequence_ready_sources": ["PDB", "NPInter"],
        "held_sources": ["RNAInter", "IntAct", "ENCODE/RBNS", "ENCODE/eCLIP (site asset pending sequence mapping)"],
        "raw_records": int(len(records)),
        "eligible_evidence_records": int(records["candidate_eligible"].sum()),
        "unique_exact_pairs": int(len(pairs)),
        "unique_proteins": int(pairs["protein_sequence"].nunique()),
        "unique_rnas": int(pairs["rna_sequence"].nunique()),
        "source_counts": {str(key): int(value) for key, value in pairs["source_category"].value_counts().items()},
        "evidence_record_count_sum": int(pairs["evidence_record_count"].sum()),
        "source_funnel": {str(row.input_source_id): {"raw_rows": int(row.raw_rows), "double_exact_rows": int(row.double_exact_rows), "eligible_evidence_rows": int(row.eligible_evidence_rows), "unique_eligible_pairs": int(row.unique_eligible_pairs)} for row in funnel.itertuples()},
        "proposal_100k_ready": bool(len(pairs) >= 100000 and not (caps["status"].eq("FAIL") | caps["status"].eq("BLOCKED")).any()),
    }
    (output / "proposal_run_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    report = [
        "# PRI-General 四主来源方案试跑报告",
        "",
        "本次运行只将两端序列和两端 exact mapping 都存在的记录纳入 pair pool；RNAInter、IntAct、ENCODE/RBNS 仅做来源审计，不补造序列。",
        "",
        "## 结果",
        "",
        f"- 原始记录：{len(records):,}",
        f"- 通过双侧 exact mapping、字母表、长度和 PDB 接触检查的证据记录：{int(records['candidate_eligible'].sum()):,}",
        f"- exact sequence pair 去重后：{len(pairs):,}",
        f"- 唯一蛋白：{pairs['protein_sequence'].nunique():,}；唯一 RNA：{pairs['rna_sequence'].nunique():,}",
        f"- 保留多条证据但每个 pair 只计一次；pair pool 的 evidence record 总数：{int(pairs['evidence_record_count'].sum()):,}",
        "",
        "## NPInter resolver 漏斗",
        "",
        markdown_table(funnel),
        "",
        "NPInter main 的解析来源组合（前 20 个）：",
        "",
        markdown_table(resolution_profile.head(20)),
        "",
        "## 来源可用性",
        "",
        markdown_table(inventory[["source_category", "candidate_pool_status", "observed_assets", "observed_rows", "local_path", "note"]]),
        "",
        "## 方案门槛",
        "",
        markdown_table(caps),
        "",
        "结论：来源比例上限本次仅作观察项，不作为硬门禁。当前新增来源虽已具备原始资产（eCLIP 也已下载），但 RNAInter/IntAct 尚未满足 exact protein + exact RNA，eCLIP 仍需基因组坐标到 RNA 窗口的解析；同时当前项目没有 Protein40/RNA80/Rfam 赋值，严格家族隔离和后续 Boltz-2 仍需补齐。",
    ]
    (output / "proposal_run_report.md").write_text("\n".join(report), encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--external-source-root", type=Path, default=None, help="local cache containing RNAInter/IntAct/legacy ENCODE assets; defaults to data/downloads")
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args()
    root = args.root.resolve()
    output = (args.output_dir or root / "reports" / "final" / "rebuild" / "proposal_run").resolve()
    records = load_candidates(root)
    pairs, evidence = make_pair_pool(records)
    stages = stage_summary(records, pairs)
    funnel = source_funnel(records)
    resolution_profile = npinter_resolution_profile(records)
    external_root = (args.external_source_root or root / "data" / "downloads").resolve()
    inventory = pd.DataFrame(source_inventory(root, external_root, records))
    caps = cap_audit(pairs)
    summary = write_report(output, records, pairs, stages, funnel, resolution_profile, inventory, caps)
    print(json.dumps({"output_dir": str(output), **summary}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
