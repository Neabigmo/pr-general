from __future__ import annotations

from collections import Counter
from pathlib import Path

from .constants import CANONICAL_FIELDS, PROTEIN_ALPHABET, RNA_ALPHABET
from .constructs import validate_construct
from .evidence import evidence_is_strict_experimental
from .families import missing_family_fields
from .split import validate_split_isolation


class ValidationError(ValueError):
    pass


def validate_records(records: list[dict], config: dict, *, strict: bool = False) -> dict:
    scope = config["scope"]
    errors: list[str] = []
    ids: set[str] = set()
    for index, record in enumerate(records):
        missing = [field for field in CANONICAL_FIELDS if field not in record]
        if missing:
            errors.append(f"row {index} missing fields: {','.join(missing)}")
            continue
        sample_id = str(record.get("sample_id") or "")
        if not sample_id or sample_id in ids:
            errors.append(f"duplicate/empty sample_id at row {index}")
        ids.add(sample_id)
        protein = str(record.get("protein_sequence") or "")
        rna = str(record.get("rna_sequence") or "")
        if set(protein) - PROTEIN_ALPHABET:
            errors.append(f"invalid protein alphabet: {sample_id}")
        if set(rna) - RNA_ALPHABET:
            errors.append(f"invalid RNA alphabet: {sample_id}")
        errors.extend(f"{sample_id}:{reason}" for reason in validate_construct(record, scope))
        if strict and record.get("split") in {"validation", "test"}:
            if not evidence_is_strict_experimental(record):
                errors.append(f"non-experimental {record['split']}: {sample_id}")
            missing_families = missing_family_fields(record)
            if missing_families:
                errors.append(f"missing holdout families {sample_id}: {','.join(missing_families)}")
    by_split = Counter(str(row.get("split") or "train") for row in records)
    target = int(scope["training_mother_samples"])
    if strict and by_split["train"] != target:
        errors.append(f"training mother sample count {by_split['train']} != {target}")
    if strict:
        try:
            validate_split_isolation(records)
        except ValueError as exc:
            errors.append(str(exc))
    if errors:
        raise ValidationError("; ".join(errors[:20]))
    return {"rows": len(records), "splits": dict(by_split), "strict": strict}


def validate_source_manifest(manifest: dict, repo_root: str | Path, *, strict: bool = False) -> list[str]:
    root = Path(repo_root)
    issues: list[str] = []
    for source in manifest.get("sources", []):
        status = source.get("status")
        local_path = source.get("local_path")
        if status == "available" and local_path and not (root / local_path).exists():
            issues.append(f"{source.get('id')}: declared available but missing {local_path}")
        if strict and status != "available" and source.get("role") != "external_structure_prediction":
            issues.append(f"{source.get('id')}: status={status}")
    return issues
