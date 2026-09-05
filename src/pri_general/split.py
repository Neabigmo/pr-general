from __future__ import annotations

import hashlib
from collections import defaultdict

from .evidence import evidence_is_strict_experimental


def _key(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def assign_splits(records: list[dict], validation_fraction: float, test_fraction: float) -> list[dict]:
    if validation_fraction < 0 or test_fraction < 0 or validation_fraction + test_fraction >= 1:
        raise ValueError("invalid split fractions")
    groups: dict[str, list[dict]] = defaultdict(list)
    for record in records:
        groups[str(record.get("biological_pair_id") or record.get("sample_id"))].append(record)
    family_groups: dict[tuple[str, str], list[str]] = defaultdict(list)
    for group_id, group in groups.items():
        experimental = next((r for r in group if evidence_is_strict_experimental(r)), None)
        if experimental is None:
            continue
        p30 = str(experimental.get("protein30") or "")
        rfam = str(experimental.get("rfam_family") or "")
        if p30 and rfam:
            family_groups[(p30, rfam)].append(group_id)
    ordered = sorted(family_groups, key=lambda pair: _key("|".join(pair)))
    test_n = round(len(ordered) * test_fraction)
    val_n = round(len(ordered) * validation_fraction)
    test_pairs = set(ordered[:test_n])
    val_pairs = set(ordered[test_n:test_n + val_n])
    assignments: dict[str, str] = {}
    for pair, group_ids in family_groups.items():
        split = "test" if pair in test_pairs else "validation" if pair in val_pairs else "train"
        assignments.update({group_id: split for group_id in group_ids})
    for group_id, group in groups.items():
        split = assignments.get(group_id, "train")
        for record in group:
            record["split"] = split
    return records


def validate_split_isolation(records: list[dict]) -> None:
    by_split: dict[str, list[dict]] = defaultdict(list)
    for record in records:
        by_split[str(record.get("split") or "train")].append(record)
    for left, right in (("train", "validation"), ("train", "test"), ("validation", "test")):
        left_pairs = {r.get("biological_pair_id") for r in by_split[left]}
        right_pairs = {r.get("biological_pair_id") for r in by_split[right]}
        if left_pairs & right_pairs:
            raise ValueError(f"biological pair leakage between {left} and {right}")
    train_p30 = {r.get("protein30") for r in by_split["train"] if r.get("protein30")}
    train_rfam = {r.get("rfam_family") for r in by_split["train"] if r.get("rfam_family")}
    for split in ("validation", "test"):
        if any(not evidence_is_strict_experimental(r) for r in by_split[split]):
            raise ValueError(f"{split} contains non-experimental evidence")
        if train_p30 & {r.get("protein30") for r in by_split[split] if r.get("protein30")}:
            raise ValueError(f"Protein30 leakage into {split}")
        if train_rfam & {r.get("rfam_family") for r in by_split[split] if r.get("rfam_family")}:
            raise ValueError(f"Rfam leakage into {split}")
