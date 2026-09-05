from __future__ import annotations

from pathlib import Path

from .io import iter_table

# Boltz-2 is executed on the school host. This module only checks the result
# shape when those externally produced results are brought back to a release.
REQUIRED_RESULT_FIELDS = ("sample_id", "confidence", "structure_path")


def validate_external_result(record: dict) -> list[str]:
    errors = [field for field in REQUIRED_RESULT_FIELDS if not record.get(field)]
    for field in ("confidence",):
        try:
            value = float(record.get(field))
        except (TypeError, ValueError):
            errors.append(f"invalid_{field}")
        else:
            if not 0 <= value <= 1:
                errors.append(f"invalid_{field}")
    return errors


def load_external_results(path: str | Path) -> list[dict]:
    rows = list(iter_table(path))
    failures = [(row.get("sample_id"), validate_external_result(row)) for row in rows]
    failures = [(sample_id, errors) for sample_id, errors in failures if errors]
    if failures:
        raise ValueError(f"invalid external Boltz result rows: {failures[:3]}")
    return rows
