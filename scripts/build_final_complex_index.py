"""Build one row per final PRI-General complex from the final audit table.

The current release counts complexes, not technical source rows:
* RCSB_PDB: one ``pdb_id`` is one experimental complex;
* other sources: one ``sequence_hash`` is one predicted complex.

The audit table is retained for provenance. This script creates a compact
index that preserves all distinct protein/RNA components for PDB complexes.
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import pandas as pd

from common import clean


ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "data" / "final" / "audit" / "general_stage_audit.csv"
OUT = ROOT / "data" / "final" / "complexes"


def true_value(value: object) -> bool:
    return clean(value).lower() in {"true", "1", "yes", "y"}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    columns = [
        "source_database", "pdb_id", "sequence_hash", "protein_sequence", "rna_sequence",
        "protein_name", "protein_accession", "rna_accession", "rna_name", "organism",
        "taxonomy_id", "system_class", "system_subtype", "rna_class", "interaction_evidence",
        "experimental_structure_available", "sample_id", "stage_queue",
    ]
    records: dict[str, dict] = {}
    for chunk in pd.read_csv(AUDIT, usecols=columns, chunksize=50_000, low_memory=False):
        chunk = chunk.loc[:, columns]
        for row in chunk.itertuples(index=False, name=None):
            (
                source_raw, pdb_id_raw, pair_hash_raw, protein_raw, rna_raw,
                protein_name_raw, protein_accession_raw, rna_accession_raw, rna_name_raw,
                organism_raw, taxid_raw, system_raw, subtype_raw, rna_class_raw,
                evidence_raw, structure_raw, sample_id_raw, queue_raw,
            ) = row
            if not true_value(queue_raw):
                continue
            source_raw = clean(source_raw)
            pdb_id = clean(pdb_id_raw)
            pair_hash = clean(pair_hash_raw)
            protein = clean(protein_raw)
            rna = clean(rna_raw)
            if not protein or not rna:
                continue
            if source_raw == "RCSB_PDB" and pdb_id:
                complex_id = f"PDB:{pdb_id}"
                complex_type = "experimental_structure"
            else:
                complex_id = f"PAIR:{pair_hash}"
                complex_type = "prediction_candidate"
            rec = records.setdefault(
                complex_id,
                {
                    "complex_id": complex_id,
                    "complex_type": complex_type,
                    "source_database": source_raw,
                    "pdb_id": pdb_id,
                    "sequence_hash": pair_hash,
                    "protein_sequences": set(),
                    "rna_sequences": set(),
                    "sample_ids": set(),
                    "protein_names": set(),
                    "protein_accessions": set(),
                    "rna_names": set(),
                    "rna_accessions": set(),
                    "organisms": set(),
                    "taxonomy_ids": set(),
                    "system_classes": set(),
                    "system_subtypes": set(),
                    "rna_classes": set(),
                    "interaction_evidence": set(),
                    "experimental_structure_available": set(),
                },
            )
            rec["protein_sequences"].add(protein)
            rec["rna_sequences"].add(rna)
            for field, value in {
                "sample_ids": sample_id_raw,
                "protein_names": protein_name_raw,
                "protein_accessions": protein_accession_raw,
                "rna_names": rna_name_raw,
                "rna_accessions": rna_accession_raw,
                "organisms": organism_raw,
                "taxonomy_ids": taxid_raw,
                "system_classes": system_raw,
                "system_subtypes": subtype_raw,
                "rna_classes": rna_class_raw,
                "interaction_evidence": evidence_raw,
                "experimental_structure_available": structure_raw,
            }.items():
                value = clean(value)
                if value:
                    rec[field].add(value)

    rows = []
    for rec in records.values():
        proteins = sorted(rec.pop("protein_sequences"))
        rnas = sorted(rec.pop("rna_sequences"))
        row = {
            "complex_id": rec["complex_id"],
            "complex_type": rec["complex_type"],
            "source_database": rec["source_database"],
            "pdb_id": rec["pdb_id"],
            "sequence_hash": rec["sequence_hash"],
            "protein_component_count": len(proteins),
            "rna_component_count": len(rnas),
            "protein_sequences_json": json.dumps(proteins, ensure_ascii=False, separators=(",", ":")),
            "rna_sequences_json": json.dumps(rnas, ensure_ascii=False, separators=(",", ":")),
        }
        for field in [
            "sample_ids", "protein_names", "protein_accessions", "rna_names", "rna_accessions",
            "organisms", "taxonomy_ids", "system_classes", "system_subtypes", "rna_classes",
            "interaction_evidence", "experimental_structure_available",
        ]:
            row[field] = ";".join(sorted(rec[field]))
        rows.append(row)

    frame = pd.DataFrame(rows).sort_values(["complex_type", "source_database", "complex_id"]).reset_index(drop=True)
    frame.insert(0, "complex_row_number", range(1, len(frame) + 1))
    frame.to_parquet(OUT / "PRI-General_final_complexes.parquet", index=False)
    frame.to_csv(OUT / "PRI-General_final_complexes.csv", index=False, encoding="utf-8-sig")
    frame[frame["complex_type"].eq("experimental_structure")].to_csv(OUT / "pdb_final_complexes.csv", index=False, encoding="utf-8-sig")
    summary = {
        "total_complexes": int(len(frame)),
        "experimental_structure_complexes": int(frame["complex_type"].eq("experimental_structure").sum()),
        "prediction_candidate_complexes": int(frame["complex_type"].eq("prediction_candidate").sum()),
        "source_counts": {str(k): int(v) for k, v in frame["source_database"].value_counts().items()},
        "statistical_unit": "PDB unique pdb_id; non-PDB unique sequence_hash",
        "input_audit": str(AUDIT),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
