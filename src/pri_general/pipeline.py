from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from typing import Any

from .config import load_config
from .constructs import validate_construct
from .ingest import load_records
from .interface import quality_tier, qualify_interface
from .io import write_json, write_jsonl
from .release import write_checksums
from .validate import validate_source_manifest

# The official NPInter main table has IDs/evidence but no sequences. It is
# recorded in provenance and resolved in a later annotation pass, not loaded
# into the sequence candidate table prematurely.
BUILD_SOURCE_IDS = ("npinter_bindingsite", "npinter_mirna", "rcsb_pdb_candidates")


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


def _annotate(record: dict, config: dict) -> dict:
    record["construct_errors"] = validate_construct(record, config["scope"])
    if record["evidence_type"] == "experimental_structure":
        qualified, reason = qualify_interface(record, config["structure"])
        record["interface_qc_reason"] = reason
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
    for source_id in BUILD_SOURCE_IDS:
        source = sources.get(source_id)
        if not source:
            skipped[source_id] = "not_in_manifest"
            continue
        path = source.get("local_path")
        if source.get("status") != "available" or not path:
            skipped[source_id] = str(source.get("status"))
            continue
        source_rows = load_records(root / path, source_id)
        rows.extend(_annotate(row, config) for row in source_rows)
        source_counts[source_id] = len(source_rows)
    skipped["npinter_main"] = "requires_sequence_resolution"

    cache_root = root / config["paths"]["build_cache"]
    candidate_path = cache_root / "candidates.jsonl"
    candidate_count = write_jsonl(candidate_path, rows)
    blockers = validate_source_manifest(source_manifest, root, strict=False)
    blockers.extend(
        f"source_{item.get('id')}_status_{item.get('status')}"
        for item in source_manifest.get("sources", [])
        if item.get("status") not in {"available", "external_host_out_of_scope"}
    )
    blockers.append("source_npinter_main_requires_sequence_resolution")
    blockers.extend([
        "family_assignments_required_before_selection",
        "external_boltz_results_required_before_final_release",
        "MSA_and_Boltz_execution_are_external_scope",
    ])
    if any(row.get("evidence_type") == "experimental_structure" and not row.get("interface_qualified") for row in rows):
        blockers.append("experimental_interface_metrics_incomplete_for_some_rows")
    stats = {
        "candidate_rows": candidate_count,
        "source_rows": dict(source_counts),
        "skipped_sources": skipped,
        "evidence_types": dict(Counter(row["evidence_type"] for row in rows)),
        "family_ready_rows": sum(bool(row["family_ready"]) for row in rows),
        "interface_qualified_rows": sum(bool(row.get("interface_qualified")) for row in rows),
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
