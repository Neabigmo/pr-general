from __future__ import annotations

import subprocess
from typing import Iterable

FAMILY_FIELDS = ("protein90", "protein40", "protein30", "rna90", "rna80", "rfam_family")


def missing_family_fields(record: dict) -> list[str]:
    return [field for field in FAMILY_FIELDS if not str(record.get(field) or "").strip()]


def family_ready(record: dict) -> bool:
    return not missing_family_fields(record)


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
