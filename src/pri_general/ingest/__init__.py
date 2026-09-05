from __future__ import annotations

from pathlib import Path

from ..io import iter_table
from ..normalize import normalize_record


def load_records(path: str | Path, source_id: str) -> list[dict]:
    return [normalize_record(raw, source_id) for raw in iter_table(path)]
