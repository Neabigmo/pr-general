from __future__ import annotations

from typing import Any


def _as_int(value: Any) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def validate_construct(record: dict, scope: dict) -> list[str]:
    errors: list[str] = []
    protein_length = _as_int(record.get("protein_length")) or 0
    rna_length = _as_int(record.get("rna_length")) or 0
    if not scope["protein_length_min"] <= protein_length <= scope["protein_length_max"]:
        errors.append("protein_length")
    if not scope["rna_length_min"] <= rna_length <= scope["rna_length_max"]:
        errors.append("rna_length")
    if record.get("site_mode") == "site_window":
        start = _as_int(record.get("binding_site_start"))
        end = _as_int(record.get("binding_site_end"))
        if start is None or end is None or end < start:
            errors.append("invalid_site_window")
        elif end - start + 1 > scope["site_window_hard_max"]:
            errors.append("site_window_too_long")
    if record.get("evidence_type") == "computational_prediction" and record.get("site_mode") != "site_window":
        errors.append("prediction_without_site_window")
    if estimate_tokens(record) > scope["token_hard_max"]:
        errors.append("token_hard_max")
    return errors


def estimate_tokens(record: dict) -> int:
    # The exact tokenizer belongs to the external Boltz host. This cheap upper
    # bound is only a preflight filter and is deliberately not a model score.
    return int(record.get("protein_length") or 0) + int(record.get("rna_length") or 0) + 16


def centered_window(length: int, center: int, target: int = 96, hard_max: int = 160) -> tuple[int, int]:
    if length <= 0 or target <= 0 or target > hard_max:
        raise ValueError("invalid sequence/window policy")
    width = min(length, target)
    start = max(1, center - (width // 2))
    end = min(length, start + width - 1)
    start = max(1, end - width + 1)
    return start, end
