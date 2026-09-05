from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[1] / "src"))

from pri_general.constants import clean_sequence
from pri_general.constructs import centered_window
from pri_general.evidence import infer_evidence_type
from pri_general.interface import qualify_interface
from pri_general.normalize import normalize_record
from pri_general.split import assign_splits, validate_split_isolation
from pri_general.validate import validate_records


def policy() -> dict:
    return {
        "scope": {
            "training_mother_samples": 1,
            "protein_length_min": 2,
            "protein_length_max": 2000,
            "rna_length_min": 2,
            "rna_length_max": 500,
            "site_window_hard_max": 160,
            "token_hard_max": 1700,
        },
        "structure": {
            "min_contact_pairs": 6,
            "min_contact_protein_residues": 3,
            "min_contact_rna_nucleotides": 3,
            "gold_resolution_max_angstrom": 3.5,
            "silver_resolution_max_angstrom": 4.0,
        },
    }


class V2CoreTests(unittest.TestCase):
    def test_sequence_normalization(self) -> None:
        self.assertEqual(clean_sequence(" a c t ", rna=True), "ACU")
        self.assertEqual(infer_evidence_type({"interaction_evidence": "E2:CLIP binding window"}), "site_resolved_binding")

    def test_centered_window_is_bounded(self) -> None:
        self.assertEqual(centered_window(200, 100, target=96), (52, 147))
        self.assertEqual(centered_window(40, 20, target=96), (1, 40))

    def test_interface_qc_does_not_guess_missing_metrics(self) -> None:
        record = {"interface_contacts": 10, "contact_protein_residues": 5, "contact_rna_nucleotides": 5}
        self.assertEqual(qualify_interface(record, policy()["structure"]), (False, "missing_backbone_atom_qc"))

    def test_split_isolates_biological_pair_and_families(self) -> None:
        records = []
        for i in range(20):
            records.append({
                "sample_id": f"s{i}", "biological_pair_id": f"bp{i}",
                "protein30": f"p{i}", "rfam_family": f"r{i}",
                "evidence_type": "experimental_structure", "split": None,
            })
        assign_splits(records, validation_fraction=0.2, test_fraction=0.2)
        validate_split_isolation(records)
        self.assertTrue({row["split"] for row in records} <= {"train", "validation", "test"})

    def test_normalized_record_is_validatable(self) -> None:
        record = normalize_record({
            "protein_sequence": "ACDE", "rna_sequence": "ACGU",
            "protein_accession": "P1", "rna_accession": "R1",
            "source_record_id": "x1", "source_database": "fixture",
            "interaction_evidence": "direct biochemical binding",
        }, "fixture")
        config = policy()
        config["scope"]["training_mother_samples"] = 1
        record.update({field: "fixture" for field in ("protein90", "protein40", "protein30", "rna90", "rna80", "rfam_family")})
        record["split"] = "train"
        result = validate_records([record], config, strict=True)
        self.assertEqual(result["rows"], 1)


if __name__ == "__main__":
    unittest.main()
