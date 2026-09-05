from __future__ import annotations

from typing import Any


def _number(value: Any) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def qualify_interface(record: dict, structure_policy: dict) -> tuple[bool, str]:
    required = {
        "interface_contacts": structure_policy["min_contact_pairs"],
        "contact_protein_residues": structure_policy["min_contact_protein_residues"],
        "contact_rna_nucleotides": structure_policy["min_contact_rna_nucleotides"],
    }
    for field, minimum in required.items():
        value = _number(record.get(field))
        if value is None:
            return False, f"missing_{field}"
        if value < minimum:
            return False, f"below_{field}"
    if record.get("backbone_atom_complete") is not True:
        return False, "missing_backbone_atom_qc"
    return True, "qualified"


def quality_tier(record: dict, structure_policy: dict) -> str | None:
    qualified, _ = qualify_interface(record, structure_policy)
    if not qualified:
        return None
    method = str(record.get("experimental_method") or "").lower()
    resolution = _number(record.get("resolution_angstrom"))
    if method in {"nmr", "solution nmr"}:
        return "gold"
    if resolution is not None and resolution <= structure_policy["gold_resolution_max_angstrom"]:
        return "gold"
    if resolution is not None and resolution <= structure_policy["silver_resolution_max_angstrom"]:
        return "silver"
    return "silver_qc"
