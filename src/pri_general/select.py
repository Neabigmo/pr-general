from __future__ import annotations

from collections import Counter

from .constructs import estimate_tokens
from .evidence import evidence_is_training_eligible, evidence_priority, evidence_support


def _quality_key(record: dict) -> tuple:
    return (
        -evidence_priority(str(record.get("evidence_type") or "")),
        -evidence_support(record),
        0 if record.get("quality_tier") == "gold" else 1,
        str(record.get("protein30") or ""),
        str(record.get("rfam_family") or ""),
        str(record.get("sample_id") or ""),
    )


def select_records(records: list[dict], target: int, caps: dict) -> tuple[list[dict], dict]:
    """Select by quality, evidence, diversity, then stable ID with hard caps."""
    selected: list[dict] = []
    counters: Counter[tuple[str, str]] = Counter()
    for record in sorted((r for r in records if evidence_is_training_eligible(r)), key=_quality_key):
        keys = _cap_keys(record)
        if any(caps.get(name) is not None and counters[(name, value)] >= int(caps[name]) for name, value in keys):
            continue
        selected.append(record)
        for key in keys:
            counters[key] += 1
        if len(selected) == target:
            break
    if len(selected) != target:
        raise ValueError(f"training target unavailable: requested {target}, selected {len(selected)}")
    return selected, {"selected": len(selected), "cap_keys": len(counters)}


def _cap_keys(record: dict) -> list[tuple[str, str]]:
    keys = [
        ("protein_exact", str(record.get("protein_sequence") or "")),
        ("protein30", str(record.get("protein30") or "")),
        ("rna_exact", str(record.get("rna_sequence") or "")),
        ("rfam", str(record.get("rfam_family") or "unknown")),
        ("source", str(record.get("source_database") or "unknown")),
        ("species", str(record.get("organism") or "unknown")),
        ("rna_class", str(record.get("rna_class") or "unknown")),
    ]
    if record.get("site_mode") == "site_window":
        keys.append(("site_window", "all"))
    organism = str(record.get("organism") or "").lower()
    if "homo sapiens" in organism or "mus musculus" in organism:
        keys.append(("human_mouse", "human_mouse"))
    if not record.get("rfam_family") or not record.get("protein30"):
        keys.append(("unknown", "unknown"))
    if estimate_tokens(record) > 1200:
        keys.append(("expensive", "expensive"))
    return [(name, value) for name, value in keys if value]
