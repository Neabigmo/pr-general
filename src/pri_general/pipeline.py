from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from typing import Any

from .config import load_config
from .constants import PROTEIN_ALPHABET, RNA_ALPHABET
from .constructs import validate_construct
from .evidence import evidence_is_strict_experimental
from .families import annotate_family_records, training_ready
from .ingest import load_records
from .interface import quality_tier, qualify_interface
from .io import write_json, write_jsonl
from .release import write_checksums
from .select import SelectionError, select_records
from .split import assign_splits, holdout_family_keys
from .structure_qc import apply_structure_qc
from .validate import validate_source_manifest

# NPInter main is loaded from its locally resolved, auditable subset; the raw
# ZIP remains the source-of-record in the manifest.
BUILD_SOURCE_IDS = ("npinter_main", "npinter_bindingsite", "npinter_mirna", "rcsb_pdb_candidates")


def _deduplicate_records(rows: list[dict]) -> tuple[list[dict], int]:
    seen: set[str] = set()
    unique: list[dict] = []
    for row in rows:
        key = json.dumps(row, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
        if key in seen:
            continue
        seen.add(key)
        unique.append(row)
    return unique, len(rows) - len(unique)


def _candidate_rejection_reasons(record: dict) -> list[str]:
    reasons = list(record.get("construct_errors") or [])
    if set(record.get("protein_sequence") or "") - PROTEIN_ALPHABET:
        reasons.append("invalid_protein_alphabet")
    if set(record.get("rna_sequence") or "") - RNA_ALPHABET:
        reasons.append("invalid_rna_alphabet")
    return sorted(set(reasons))


def load_source_manifest(path: str | Path) -> dict[str, Any]:
    try:
        import yaml
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("PyYAML is required to load source provenance") from exc
    data = yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}
    if data.get("schema_version") != "pri-general/source-manifest/v2":
        raise ValueError("invalid source manifest schema")
    return data


def _source_map(manifest: dict) -> dict[str, dict]:
    return {str(item["id"]): item for item in manifest.get("sources", [])}


def _source_inventory(manifest: dict) -> dict[str, dict]:
    return {
        str(item["id"]): {
            "role": item.get("role"),
            "status": item.get("status"),
            "local_path": item.get("local_path"),
            "resolved_path": item.get("resolved_path"),
            "disposition": item.get("disposition", "provenance_only"),
        }
        for item in manifest.get("sources", [])
    }


def _annotate(record: dict, config: dict) -> dict:
    record["construct_errors"] = validate_construct(record, config["scope"])
    if record["evidence_type"] == "experimental_structure":
        prior_qc_reason = record.get("interface_qc_reason")
        qualified, reason = qualify_interface(record, config["structure"])
        record["interface_qc_reason"] = (
            prior_qc_reason
            if prior_qc_reason and prior_qc_reason != "qualified_metrics"
            else reason
        )
        record["quality_tier"] = quality_tier(record, config["structure"])
        record["interface_qualified"] = qualified
    else:
        record["interface_qc_reason"] = None
        record["interface_qualified"] = None
    record["family_ready"] = not any(not record.get(field) for field in ("protein30", "rna80", "rfam_family"))
    return record


def build(repo_root: str | Path, config_path: str | Path, release: str) -> dict:
    root = Path(repo_root)
    config = load_config(config_path)
    manifest_path = root / config["paths"]["source_manifest"]
    source_manifest = load_source_manifest(manifest_path)
    sources = _source_map(source_manifest)
    rows: list[dict] = []
    source_counts: Counter[str] = Counter()
    skipped: dict[str, str] = {}
    cache_root = root / config["paths"]["build_cache"]
    for source_id in BUILD_SOURCE_IDS:
        source = sources.get(source_id)
        if not source:
            skipped[source_id] = "not_in_manifest"
            continue
        path = source.get("local_path")
        if source.get("status") != "available" or not path:
            skipped[source_id] = str(source.get("status"))
            continue
        resolved_path = source.get("resolved_path")
        input_path = root / (resolved_path or path)
        if not input_path.exists():
            skipped[source_id] = f"missing_resolved_input:{input_path}"
            continue
        source_rows = load_records(input_path, source_id)
        rows.extend(source_rows)
        source_counts[source_id] = len(source_rows)
    rows, duplicate_rows_removed = _deduplicate_records(rows)
    structure_qc_count = apply_structure_qc(rows, cache_root / "structure_qc.jsonl")
    rows = [_annotate(row, config) for row in rows]
    rejected_rows: list[dict] = []
    accepted_rows: list[dict] = []
    rejection_reasons: Counter[str] = Counter()
    for row in rows:
        reasons = _candidate_rejection_reasons(row)
        if reasons:
            rejected_rows.append(row)
            rejection_reasons.update(reasons)
        else:
            accepted_rows.append(row)
    rejection_sources = Counter(row["source_id"] for row in rejected_rows)
    rows = accepted_rows
    skipped.update({
        source_id: source.get("disposition", "provenance_only")
        for source_id, source in sources.items()
        if source_id not in BUILD_SOURCE_IDS
    })
    family_stats = annotate_family_records(rows, cache_root / "family_work")
    candidate_path = cache_root / "candidates.jsonl"
    holdout_source = [
        row for row in rows
        if evidence_is_strict_experimental(row) and row.get("interface_qualified") and training_ready(row)
    ]
    assign_splits(
        holdout_source,
        validation_fraction=float(config["splits"]["validation_fraction"]),
        test_fraction=float(config["splits"]["test_fraction"]),
    )
    holdout_pairs = {
        row["biological_pair_id"] for row in holdout_source if row.get("split") in {"validation", "test"}
    }
    holdout_keys = {
        key
        for row in holdout_source
        if row.get("split") in {"validation", "test"}
        for key in holdout_family_keys(row)
    }
    train_pool = [
        row for row in rows
        if (
            training_ready(row)
            and row["biological_pair_id"] not in holdout_pairs
            and not holdout_family_keys(row) & holdout_keys
        )
    ]
    selected_train: list[dict] = []
    selection_error = None
    try:
        selected_train, selection_stats = select_records(
            train_pool,
            target=int(config["scope"]["training_mother_samples"]),
            caps=config["caps"],
        )
        for row in selected_train:
            row["split"] = "train"
    except SelectionError as exc:
        selected_train = exc.selected
        for row in selected_train:
            row["split"] = "train"
        selection_stats = {"selected": len(selected_train), "error": str(exc)}
        selection_error = "training_target_unavailable"
    write_jsonl(cache_root / "selected_train.jsonl", selected_train)
    write_jsonl(cache_root / "validation_candidates.jsonl", [row for row in holdout_source if row.get("split") == "validation"])
    write_jsonl(cache_root / "test_candidates.jsonl", [row for row in holdout_source if row.get("split") == "test"])
    candidate_count = write_jsonl(candidate_path, rows)
    blockers = validate_source_manifest(source_manifest, root, strict=False)
    blockers.extend(
        f"source_{item.get('id')}_status_{item.get('status')}"
        for item in source_manifest.get("sources", [])
        if item.get("status") not in {"available", "external_host_out_of_scope"}
    )
    npinter = sources.get("npinter_main", {})
    if npinter.get("disposition") == "canonical_candidate_input_resolved_subset":
        blockers.append("source_npinter_main_partial_sequence_resolution")
    else:
        blockers.append("source_npinter_main_requires_sequence_resolution")
    if family_stats["assignment_files_present"] < 6 or not family_stats["training_ready_rows"]:
        blockers.append("family_assignments_required_before_selection")
    if selection_error:
        blockers.append(selection_error)
    if not holdout_source:
        blockers.append("strict_experimental_holdout_unavailable")
    blockers.extend(["external_boltz_results_required_before_final_release", "MSA_and_Boltz_execution_are_external_scope"])
    if any(row.get("evidence_type") == "experimental_structure" and not row.get("interface_qualified") for row in rows):
        blockers.append("experimental_interface_metrics_incomplete_for_some_rows")
    stats = {
        "candidate_rows": candidate_count,
        "duplicate_rows_removed": duplicate_rows_removed,
        "rejected_rows": len(rejected_rows),
        "rejection_reasons": dict(rejection_reasons),
        "rejection_sources": dict(rejection_sources),
        "source_rows": dict(source_counts),
        "skipped_sources": skipped,
        "evidence_types": dict(Counter(row["evidence_type"] for row in rows)),
        "family_ready_rows": int(family_stats["family_ready_rows"]),
        "interface_qualified_rows": sum(bool(row.get("interface_qualified")) for row in rows),
        "structure_qc_rows_applied": structure_qc_count,
        "family_assignments": family_stats,
        "selection": selection_stats,
        "validation_candidates": sum(row.get("split") == "validation" for row in holdout_source),
        "test_candidates": sum(row.get("split") == "test" for row in holdout_source),
        "strict_holdout_source_rows": len(holdout_source),
        "source_inventory": _source_inventory(source_manifest),
    }
    release_root = root / config["paths"]["release_root"] / release
    release_root.mkdir(parents=True, exist_ok=True)
    write_json(release_root / "statistics.json", stats)
    write_json(release_root / "manifest.json", {
        "schema_version": "pri-general/release-manifest/v2",
        "release": release,
        "status": "DRAFT",
        "candidate_path": candidate_path.relative_to(root).as_posix(),
        "blockers": sorted(set(blockers)),
        "source_manifest": manifest_path.relative_to(root).as_posix(),
    })
    release_files = [
        path for path in release_root.rglob("*")
        if path.is_file() and path.name != "checksums.sha256"
    ]
    write_checksums(release_files, release_root / "checksums.sha256", base_dir=root)
    return stats | {"release_status": "DRAFT", "blockers": sorted(set(blockers))}


def report(repo_root: str | Path, config_path: str | Path, release: str) -> dict:
    root = Path(repo_root)
    config = load_config(config_path)
    release_root = root / config["paths"]["release_root"] / release
    return json.loads((release_root / "manifest.json").read_text(encoding="utf-8"))
