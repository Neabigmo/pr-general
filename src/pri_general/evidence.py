from __future__ import annotations

from .constants import EVIDENCE_TYPES

EVIDENCE_PRIORITY = {
    "experimental_structure": 5,
    "direct_biochemical": 4,
    "site_resolved_binding": 3,
    "physical_association": 2,
    "computational_prediction": 1,
}


def infer_evidence_type(raw: dict) -> str:
    text = " ".join(
        str(raw.get(key) or "").lower()
        for key in ("interaction_evidence", "evidence_type", "evidence", "experiment", "class", "source_database")
    )
    if raw.get("experimental_structure_available") or "pdb" in text or "structure" in text:
        return "experimental_structure"
    if any(token in text for token in ("clip", "binding site", "bindingsite", "peak", "site-resolved")):
        return "site_resolved_binding"
    if any(token in text for token in ("binding", "kd", "affinity", "biochemical", "native pairing")):
        return "direct_biochemical"
    if any(token in text for token in ("association", "interact", "physical")):
        return "physical_association"
    return "computational_prediction"


def evidence_priority(value: str) -> int:
    return EVIDENCE_PRIORITY.get(value, 0)


def evidence_is_training_eligible(record: dict) -> bool:
    return record.get("evidence_type") in EVIDENCE_TYPES[:-1]


def evidence_is_strict_experimental(record: dict) -> bool:
    return record.get("evidence_type") == "experimental_structure"


def evidence_support(record: dict) -> float:
    try:
        return float(record.get("evidence_support"))
    except (TypeError, ValueError):
        return float(evidence_priority(str(record.get("evidence_type") or "")))
