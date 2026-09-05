"""Portable entry point for the current PRI-General release only."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

from scripts.common import load_pipeline_config, pair_hash


ROOT = Path(__file__).resolve().parent
FINAL = ROOT / "data" / "final"


def rows(path: Path) -> int:
    return pq.ParquetFile(path).metadata.num_rows


def validate() -> dict:
    config = load_pipeline_config(ROOT)
    protein_min, protein_max = map(int, config["protein_length"])
    rna_min, rna_max = map(int, config["rna_length"])
    expected_counts = {
        "index": int(config["final_complex_count"]),
        "queue": int(config["pre_boltz_queue_count"]),
        "d1": int(config["d1_hard_qc_count"]),
        "d2": int(config["d2_exact_dedup_count"]),
        "d3": int(config["d3_local_rna_representative_count"]),
        "pdb": int(config["experimental_pdb_complex_count"]),
    }
    checks = []
    index = FINAL / "complexes" / "PRI-General_final_complexes.parquet"
    summary = FINAL / "complexes" / "summary.json"
    queue = FINAL / "pre_boltz" / "prediction_queue.parquet"
    d1 = FINAL / "pre_boltz" / "general_d1_hard_qc.parquet"
    d2 = FINAL / "pre_boltz" / "general_d2_exact_dedup.parquet"
    d3 = FINAL / "pre_boltz" / "general_d3_local_rna_representatives.parquet"
    audit = FINAL / "audit" / "general_stage_audit.csv"

    expected = {
        index: expected_counts["index"],
        queue: expected_counts["queue"],
        d1: expected_counts["d1"],
        d2: expected_counts["d2"],
        d3: expected_counts["d3"],
    }
    for path, expected_rows in expected.items():
        observed = rows(path) if path.exists() else None
        checks.append({"check": f"rows:{path.relative_to(ROOT)}", "ok": observed == expected_rows, "observed": observed, "expected": expected_rows})
    checks.append({"check": "audit_exists", "ok": audit.exists(), "observed": str(audit)})

    if index.exists():
        frame = pd.read_parquet(index, columns=["complex_id", "complex_type", "source_database"])
        checks.append({"check": "complex_id_unique", "ok": frame["complex_id"].is_unique, "observed": int(frame["complex_id"].nunique()), "expected": len(frame)})
        pdb_count = int((frame["complex_type"] == "experimental_structure").sum())
        checks.append({"check": "pdb_complex_count", "ok": pdb_count == expected_counts["pdb"], "observed": pdb_count, "expected": expected_counts["pdb"]})

    if queue.exists():
        q = pd.read_parquet(queue, columns=["sample_id", "sequence_hash", "protein_sequence", "rna_sequence", "benchmark_conflict"])
        hash_mismatch = sum(pair_hash(str(p), str(r)) != str(h) for p, r, h in zip(q["protein_sequence"], q["rna_sequence"], q["sequence_hash"]))
        invalid = ((q["protein_sequence"].astype(str).str.len() < protein_min) | (q["protein_sequence"].astype(str).str.len() > protein_max) | (q["rna_sequence"].astype(str).str.len() < rna_min) | (q["rna_sequence"].astype(str).str.len() > rna_max) | q["benchmark_conflict"].fillna(False).astype(bool)).sum()
        checks.extend([
            {"check": "queue_sample_id_unique", "ok": q["sample_id"].is_unique, "observed": int(q["sample_id"].nunique()), "expected": len(q)},
            {"check": "queue_sequence_hash_unique", "ok": q["sequence_hash"].is_unique, "observed": int(q["sequence_hash"].nunique()), "expected": len(q)},
            {"check": "queue_hash_matches_sequence", "ok": hash_mismatch == 0, "observed": hash_mismatch, "expected": 0},
            {"check": "queue_length_and_benchmark_qc", "ok": int(invalid) == 0, "observed": int(invalid), "expected": 0},
        ])

    result = {"dataset": "PRI-General", "root": str(ROOT), "passed": all(item["ok"] for item in checks), "checks": checks}
    report = ROOT / "reports" / "final" / "validation.json"
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result


def run_script(name: str) -> None:
    subprocess.run([sys.executable, str(ROOT / "scripts" / name)], cwd=ROOT, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="当前 PRI-General 独立工作区入口")
    parser.add_argument("command", choices=["validate", "analyze", "rebuild", "build-index"], help="执行验证、分析、从源候选重跑或重建复合物索引")
    args = parser.parse_args()
    if args.command == "validate":
        result = validate()
        raise SystemExit(0 if result["passed"] else 1)
    if args.command == "build-index":
        run_script("build_final_complex_index.py")
    elif args.command == "analyze":
        run_script("build_final_complex_index.py")
        run_script("make_source_stage_overview_table.py")
        run_script("analyze_final_complexes.py")
        run_script("make_compact_analysis_dashboard.py")
        run_script("make_extended_analysis_figures.py")
    else:
        run_script("pri_general_full_pipeline.py")


if __name__ == "__main__":
    main()
