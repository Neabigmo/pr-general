from __future__ import annotations

import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Iterable

FAMILY_FIELDS = ("protein90", "protein40", "protein30", "rna90", "rna80", "rfam_family")


def missing_family_fields(record: dict) -> list[str]:
    return [field for field in FAMILY_FIELDS if not str(record.get(field) or "").strip()]


def family_ready(record: dict) -> bool:
    return not missing_family_fields(record)


def training_ready(record: dict) -> bool:
    """Return whether the record has the family keys needed for training.

    Rfam is the preferred RNA family key, but the workflow explicitly falls
    back to RNA80 when Rfam has no annotation.  Strict holdout validation can
    apply the same fallback without inventing a family assignment.
    """
    return all(str(record.get(field) or "").strip() for field in ("protein30", "rna80"))


def mmseqs2_version() -> str:
    completed = subprocess.run(
        ["wsl.exe", "-d", "Ubuntu", "--", "mmseqs", "version"],
        check=True, capture_output=True, text=True,
    )
    return (completed.stdout or completed.stderr).strip()


def require_family_fields(records: Iterable[dict]) -> None:
    for record in records:
        missing = missing_family_fields(record)
        if missing:
            raise ValueError(f"family assignments missing for {record.get('sample_id')}: {', '.join(missing)}")


def _cluster_map(path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    if not path.exists():
        return mapping
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[0] and parts[1]:
                mapping[parts[1]] = parts[0]
    return mapping


def _rfam_map(path: Path) -> dict[str, str]:
    values: dict[str, set[str]] = defaultdict(set)
    if not path.exists():
        return {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            raw_id, family = line.rstrip("\n").split("\t", 1)
            if raw_id and family:
                values[raw_id.split(".", 1)[0]].add(family)
    return {key: next(iter(items)) for key, items in values.items() if len(items) == 1}


def annotate_family_records(records: list[dict], work_root: str | Path) -> dict[str, int | bool]:
    """Apply local MMseqs2/RNAcentral/Rfam assignments without guessing."""
    root = Path(work_root)
    maps = {
        field: _cluster_map(root / f"{field}_cluster.tsv")
        for field in ("protein90", "protein40", "protein30", "rna90", "rna80")
    }
    rfam = _rfam_map(root / "rna_rfam_by_id.tsv")
    rfam_by_sequence = _rfam_map(root / "rna_rfam_by_seq.tsv")
    assigned = 0
    training_ready_rows = 0
    for record in records:
        sequence_key = str(record.get("sequence_key") or "")
        protein_key, separator, rna_digest = sequence_key.partition("_r_")
        rna_key = f"r_{rna_digest}" if separator else ""
        for field, mapping in maps.items():
            key = protein_key if field.startswith("protein") else rna_key
            if mapping and key:
                record[field] = mapping.get(key)
        rna_id = str(record.get("rna_id") or "").split(".", 1)[0]
        if rna_id and rna_id in rfam:
            record["rfam_family"] = rfam[rna_id]
        elif rna_key in rfam_by_sequence:
            record["rfam_family"] = rfam_by_sequence[rna_key]
        record["family_ready"] = family_ready(record)
        training_ready_rows += int(training_ready(record))
        assigned += int(record["family_ready"])
    return {
        "family_ready_rows": assigned,
        "training_ready_rows": training_ready_rows,
        "protein90_assignments": len(maps["protein90"]),
        "protein40_assignments": len(maps["protein40"]),
        "protein30_assignments": len(maps["protein30"]),
        "rna90_assignments": len(maps["rna90"]),
        "rna80_assignments": len(maps["rna80"]),
        "rfam_assignments": len(rfam),
        "rfam_sequence_assignments": len(rfam_by_sequence),
        "assignment_files_present": sum(bool(mapping) for mapping in maps.values())
        + int(bool(rfam))
        + int(bool(rfam_by_sequence)),
    }
