"""Measure sequence-resolver yield on local RNAInter and IntAct samples.

The test is intentionally bounded.  It uses accession/xref identifiers and
the sequence maps already present in the current candidate tables; names are
never treated as exact sequence mappings.
"""

from __future__ import annotations

import argparse
import json
import re
import tarfile
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd


RNA_CATEGORY_RE = re.compile(r"rna|lncrna|mirna|mrna|ncrna|snrna|snorna|rrna|trna|circrna", re.I)
PROTEIN_CATEGORY_RE = re.compile(r"protein|rbp|transcription factor|^tf$", re.I)
RNA_TYPE_RE = re.compile(r"rna|nucleic acid|mi:0327|mi:0328", re.I)
PROTEIN_TYPE_RE = re.compile(r"protein|mi:0326", re.I)


def clean(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    return "" if text.lower() in {"", "nan", "none", "null", "n/a", "na"} else text


def normalize_id(value: object) -> str:
    text = clean(value).upper()
    if not text:
        return ""
    if "|" in text:
        text = text.split("|", 1)[0]
    for prefix in ("UNIPROTKB:", "UNIPROT:", "REFSEQ:", "NCBI:", "ENSEMBL:", "MIRBASE:", "RNACENTRAL:"):
        if text.startswith(prefix):
            return text[len(prefix):]
    return text


def add_mapping(index: dict[str, set[str]], identifier: object, sequence: object, *, rna: bool = False) -> None:
    key = normalize_id(identifier)
    seq = clean(sequence).upper()
    if rna:
        seq = seq.replace("T", "U")
    if key and seq:
        index.setdefault(key, set()).add(seq)


def load_sequence_maps(root: Path) -> tuple[dict[str, set[str]], dict[str, set[str]]]:
    proteins: dict[str, set[str]] = {}
    rnas: dict[str, set[str]] = {}
    for filename in ("npinter_raw.parquet", "npinter_bindingsite_raw.parquet", "mirna_protein_raw.parquet", "pdb_raw.parquet"):
        path = root / "data" / "source_candidates" / filename
        frame = pd.read_parquet(path, columns=["protein_accession", "rna_accession", "protein_sequence", "rna_sequence"])
        for protein_id, rna_id, protein_sequence, rna_sequence in frame.itertuples(index=False, name=None):
            add_mapping(proteins, protein_id, protein_sequence)
            add_mapping(rnas, rna_id, rna_sequence, rna=True)
    return proteins, rnas


def resolve(index: dict[str, set[str]], identifier: object) -> str:
    sequences = index.get(normalize_id(identifier), set())
    return next(iter(sequences)) if len(sequences) == 1 else ""


def is_experimental(row: pd.Series) -> bool:
    strong = clean(row.get("strong")).lower()
    weak = clean(row.get("weak")).lower()
    predict = clean(row.get("predict")).lower()
    return bool((strong or weak) and not predict)


def rna_inter_sample(path: Path, sample_rows: int, proteins: dict[str, set[str]], rnas: dict[str, set[str]]) -> dict[str, object]:
    with tarfile.open(path, mode="r:gz") as archive:
        member = archive.next()
        if member is None:
            return {"archive": path.name, "raw_rows": 0, "note": "empty archive"}
        handle = archive.extractfile(member)
        if handle is None:
            return {"archive": path.name, "raw_rows": 0, "note": "member unreadable"}
        frame = pd.read_csv(handle, sep="\t", nrows=sample_rows, dtype=str, keep_default_na=False)
    category_a = frame["Category1"].astype(str)
    category_b = frame["Category2"].astype(str)
    a_is_rna = category_a.str.contains(RNA_CATEGORY_RE, na=False)
    b_is_rna = category_b.str.contains(RNA_CATEGORY_RE, na=False)
    a_is_protein = category_a.str.contains(PROTEIN_CATEGORY_RE, na=False)
    b_is_protein = category_b.str.contains(PROTEIN_CATEGORY_RE, na=False)
    a_rna_b_protein = a_is_rna & b_is_protein
    b_rna_a_protein = b_is_rna & a_is_protein
    is_pair = a_rna_b_protein | b_rna_a_protein
    experimental = frame.apply(is_experimental, axis=1)
    pair = frame.loc[is_pair].copy()
    pair["protein_id_for_resolution"] = np.where(a_rna_b_protein[is_pair], pair["Raw_ID2"], pair["Raw_ID1"])
    pair["rna_id_for_resolution"] = np.where(a_rna_b_protein[is_pair], pair["Raw_ID1"], pair["Raw_ID2"])
    pair["protein_sequence_resolved"] = pair["protein_id_for_resolution"].map(lambda value: bool(resolve(proteins, value)))
    pair["rna_sequence_resolved"] = pair["rna_id_for_resolution"].map(lambda value: bool(resolve(rnas, value)))
    pair["experimental"] = experimental[is_pair].to_numpy()
    exact = pair["experimental"] & pair["protein_sequence_resolved"] & pair["rna_sequence_resolved"]
    exact_pairs = set(zip(pair.loc[exact, "protein_id_for_resolution"], pair.loc[exact, "rna_id_for_resolution"]))
    return {
        "archive": path.name,
        "sample_rows": int(sample_rows),
        "raw_rows": int(len(frame)),
        "rna_protein_rows": int(is_pair.sum()),
        "experimental_rows": int((is_pair & experimental).sum()),
        "protein_id_present_rows": int(pair["protein_id_for_resolution"].map(clean).ne("").sum()),
        "rna_id_present_rows": int(pair["rna_id_for_resolution"].map(clean).ne("").sum()),
        "protein_sequence_resolved_rows": int((pair["experimental"] & pair["protein_sequence_resolved"]).sum()),
        "rna_sequence_resolved_rows": int((pair["experimental"] & pair["rna_sequence_resolved"]).sum()),
        "double_sequence_rows": int(exact.sum()),
        "unique_exact_identifier_pairs": len(exact_pairs),
        "note": "exact pair count is a lower-bound resolver test against local candidate-table ID maps",
    }


def intact_sample(path: Path, sample_rows: int, proteins: dict[str, set[str]], rnas: dict[str, set[str]]) -> dict[str, object]:
    with zipfile.ZipFile(path) as archive:
        with archive.open("intact.txt") as binary:
            frame = pd.read_csv(binary, sep="\t", nrows=sample_rows, dtype=str, keep_default_na=False, on_bad_lines="skip")
    frame.columns = [str(column).lstrip("#").strip() for column in frame.columns]
    type_a = frame["Type(s) interactor A"].astype(str)
    type_b = frame["Type(s) interactor B"].astype(str)
    a_is_rna = type_a.str.contains(RNA_TYPE_RE, na=False)
    b_is_rna = type_b.str.contains(RNA_TYPE_RE, na=False)
    a_is_protein = type_a.str.contains(PROTEIN_TYPE_RE, na=False)
    b_is_protein = type_b.str.contains(PROTEIN_TYPE_RE, na=False)
    a_rna_b_protein = a_is_rna & b_is_protein
    b_rna_a_protein = b_is_rna & a_is_protein
    is_pair = a_rna_b_protein | b_rna_a_protein
    interaction_type = frame["Interaction type(s)"].astype(str)
    direct = interaction_type.str.contains("MI:0407|direct interaction", case=False, regex=True, na=False)
    physical = interaction_type.str.contains("MI:0915|physical association", case=False, regex=True, na=False)
    negative = frame["Negative"].astype(str).str.lower().isin({"true", "yes", "1"})
    experimental_direct = is_pair & direct & ~physical & ~negative
    pair = frame.loc[experimental_direct].copy()
    pair["protein_id_for_resolution"] = np.where(a_rna_b_protein[experimental_direct], pair["ID(s) interactor B"], pair["ID(s) interactor A"])
    pair["rna_id_for_resolution"] = np.where(a_rna_b_protein[experimental_direct], pair["ID(s) interactor A"], pair["ID(s) interactor B"])
    pair["protein_sequence_resolved"] = pair["protein_id_for_resolution"].map(lambda value: bool(resolve(proteins, value)))
    pair["rna_sequence_resolved"] = pair["rna_id_for_resolution"].map(lambda value: bool(resolve(rnas, value)))
    exact = pair["protein_sequence_resolved"] & pair["rna_sequence_resolved"]
    exact_pairs = set(zip(pair.loc[exact, "protein_id_for_resolution"], pair.loc[exact, "rna_id_for_resolution"]))
    return {
        "archive": path.name,
        "sample_rows": int(sample_rows),
        "raw_rows": int(len(frame)),
        "protein_rna_rows": int(is_pair.sum()),
        "direct_nonphysical_rows": int(experimental_direct.sum()),
        "protein_id_present_rows": int(pair["protein_id_for_resolution"].map(clean).ne("").sum()),
        "rna_id_present_rows": int(pair["rna_id_for_resolution"].map(clean).ne("").sum()),
        "protein_sequence_resolved_rows": int((pair["protein_sequence_resolved"]).sum()),
        "rna_sequence_resolved_rows": int((pair["rna_sequence_resolved"]).sum()),
        "double_sequence_rows": int(exact.sum()),
        "unique_exact_identifier_pairs": len(exact_pairs),
        "note": "exact pair count is a lower-bound resolver test against local candidate-table ID maps",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--external-source-root", type=Path, default=None, help="local cache containing RNAInter and IntAct; defaults to data/downloads")
    parser.add_argument("--sample-rows", type=int, default=100000, help="rows sampled from each archive")
    parser.add_argument("--output-dir", type=Path, default=Path("reports/final/rebuild/external_resolver_yield"))
    args = parser.parse_args()
    root = args.root.resolve()
    external = (args.external_source_root or root / "data" / "downloads").resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    proteins, rnas = load_sequence_maps(root)
    rnainter_results = [rna_inter_sample(path, args.sample_rows, proteins, rnas) for path in sorted((external / "rnainter").glob("*.tar.gz"))]
    intact_path = external / "intact" / "intact.zip"
    intact_result = intact_sample(intact_path, args.sample_rows, proteins, rnas) if intact_path.exists() else {"archive": "intact.zip", "note": "missing"}
    result = {
        "sample_rows_per_archive": int(args.sample_rows),
        "local_protein_identifier_keys": len(proteins),
        "local_rna_identifier_keys": len(rnas),
        "rnainter": rnainter_results,
        "intact": intact_result,
        "interpretation": "A row enters a resolver yield only when the filtered experimental relation has unique local protein and RNA sequence mappings; this pass does not use gene names or create sequence windows.",
    }
    (output / "external_resolver_yield.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    rows = [{"source": "RNAInter", **row} for row in rnainter_results] + [{"source": "IntAct", **intact_result}]
    pd.DataFrame(rows).to_csv(output / "external_resolver_yield.csv", index=False, encoding="utf-8-sig")
    print(json.dumps({"output_dir": str(output), **result}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
