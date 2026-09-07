from __future__ import annotations

import hashlib
from collections import defaultdict

from .evidence import evidence_is_strict_experimental


def _key(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _rna_holdout_keys(record: dict) -> set[str]:
    keys: set[str] = set()
    rfam = str(record.get("rfam_family") or "").strip()
    if rfam:
        keys.add(f"rfam:{rfam}")
    rna80 = str(record.get("rna80") or "").strip()
    if rna80:
        keys.add(f"rna80:{rna80}")
    return keys


def holdout_family_keys(record: dict) -> set[str]:
    keys = _rna_holdout_keys(record)
    p30 = str(record.get("protein30") or "").strip()
    if p30:
        keys.add(f"protein30:{p30}")
    return keys


def assign_splits(records: list[dict], validation_fraction: float, test_fraction: float) -> list[dict]:
    if validation_fraction < 0 or test_fraction < 0 or validation_fraction + test_fraction >= 1:
        raise ValueError("invalid split fractions")
    groups: dict[str, list[dict]] = defaultdict(list)
    for record in records:
        groups[str(record.get("biological_pair_id") or record.get("sample_id"))].append(record)
    group_keys: dict[str, set[str]] = {}
    for group_id, group in groups.items():
        experimental = next((r for r in group if evidence_is_strict_experimental(r)), None)
        if experimental is None:
            continue
        p30 = str(experimental.get("protein30") or "")
        rna_keys = _rna_holdout_keys(experimental)
        if p30 and rna_keys:
            group_keys[group_id] = holdout_family_keys(experimental)

    parent = {group_id: group_id for group_id in group_keys}

    def find(group_id: str) -> str:
        while parent[group_id] != group_id:
            parent[group_id] = parent[parent[group_id]]
            group_id = parent[group_id]
        return group_id

    nodes: dict[str, str] = {}
    for group_id, keys in group_keys.items():
        for node in keys:
            previous = nodes.get(node)
            if previous is None:
                nodes[node] = group_id
                continue
            left, right = find(group_id), find(previous)
            if left != right:
                parent[right] = left

    components: dict[str, list[str]] = defaultdict(list)
    for group_id in group_keys:
        components[find(group_id)].append(group_id)
    ordered = sorted(
        components.values(),
        key=lambda group_ids: _key("|".join(sorted(group_ids))),
    )
    test_n = round(len(ordered) * test_fraction)
    val_n = round(len(ordered) * validation_fraction)
    assignments: dict[str, str] = {}
    for index, group_ids in enumerate(ordered):
        split = "test" if index < test_n else "validation" if index < test_n + val_n else "train"
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
    split_keys = {
        split: {
            "protein30": {r.get("protein30") for r in by_split[split] if r.get("protein30")},
            "rna_family": {
                key
                for r in by_split[split]
                for key in _rna_holdout_keys(r)
            },
        }
        for split in ("train", "validation", "test")
    }
    for left, right in (("train", "validation"), ("train", "test"), ("validation", "test")):
        if split_keys[left]["protein30"] & split_keys[right]["protein30"]:
            raise ValueError(f"Protein30 leakage between {left} and {right}")
        if split_keys[left]["rna_family"] & split_keys[right]["rna_family"]:
            raise ValueError(f"RNA family leakage between {left} and {right}")
    for split in ("validation", "test"):
        if any(not evidence_is_strict_experimental(r) for r in by_split[split]):
            raise ValueError(f"{split} contains non-experimental evidence")
