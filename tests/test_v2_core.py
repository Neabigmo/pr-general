from __future__ import annotations

import sys
import tempfile
from pathlib import Path
import unittest

sys.path.insert(0, str(Path(__file__).parents[1] / "src"))

from pri_general.constants import clean_sequence
from pri_general.constructs import centered_window
from pri_general.evidence import infer_evidence_type
from pri_general.interface import qualify_interface
from pri_general.normalize import normalize_record
from pri_general.pipeline import _candidate_rejection_reasons
from pri_general.families import annotate_family_records, training_ready
from pri_general.select import select_records
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

    def test_split_uses_rna80_when_rfam_is_missing(self) -> None:
        records = []
        for i in range(20):
            records.append({
                "sample_id": f"fallback{i}", "biological_pair_id": f"fallback_bp{i}",
                "protein30": f"p{i}", "rna80": f"r{i}",
                "evidence_type": "experimental_structure", "split": None,
            })
        assign_splits(records, validation_fraction=0.2, test_fraction=0.2)
        validate_split_isolation(records)
        self.assertTrue(all(training_ready(row) for row in records))

    def test_split_keeps_connected_family_components_together(self) -> None:
        records = [
            {
                "sample_id": "connected0", "biological_pair_id": "connected_bp0",
                "protein30": "p0", "rna80": "r0", "rfam_family": "RF0",
                "evidence_type": "experimental_structure", "split": None,
            },
            {
                "sample_id": "connected1", "biological_pair_id": "connected_bp1",
                "protein30": "p0", "rna80": "r1", "rfam_family": "RF1",
                "evidence_type": "experimental_structure", "split": None,
            },
        ]
        assign_splits(records, validation_fraction=0.5, test_fraction=0.0)
        self.assertEqual(records[0]["split"], records[1]["split"])
        validate_split_isolation(records)

    def test_missing_rfam_does_not_consume_rfam_cap(self) -> None:
        records = []
        for i in range(2):
            records.append({
                "sample_id": f"no_rfam{i}", "protein_sequence": "ACDE" if i == 0 else "ACDF",
                "protein30": f"p{i}", "rna_sequence": "ACGU" if i == 0 else "ACGA",
                "source_database": "fixture", "rna_class": "other_ncrna",
                "evidence_type": "direct_biochemical",
            })
        selected, _ = select_records(records, target=2, caps={"rfam": 1})
        self.assertEqual(len(selected), 2)

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

    def test_site_windows_have_distinct_sample_ids(self) -> None:
        common = {
            "protein_sequence": "ACDE",
            "rna_sequence": "ACGU",
            "protein_accession": "P1",
            "rna_accession": "R1",
            "source_record_id": "x1",
            "binding_site_start": 10,
            "source_database": "fixture",
            "interaction_evidence": "site binding window",
        }
        first = normalize_record({**common, "binding_site_end": 20}, "fixture")
        second = normalize_record({**common, "binding_site_end": 21}, "fixture")
        self.assertNotEqual(first["sample_id"], second["sample_id"])

    def test_invalid_candidate_is_rejected_before_release(self) -> None:
        record = {"construct_errors": [], "protein_sequence": "ACDE", "rna_sequence": "IIII"}
        self.assertEqual(_candidate_rejection_reasons(record), ["invalid_rna_alphabet"])

    def test_structure_notes_are_parsed_without_inventing_interface_qc(self) -> None:
        record = normalize_record({
            "protein_sequence": "ACDE", "rna_sequence": "ACGU",
            "protein_accession": "1ABC:E", "rna_accession": "1ABC:A",
            "source_record_id": "x1", "pdb_id": "1ABC",
            "experimental_structure_available": True,
            "interaction_evidence": "E3:PDB experimental complex, 12 contact atom pairs <=5A",
            "notes": "protein_chain=E; rna_chain=A; resolution=[2.38]; method=X",
        }, "fixture")
        self.assertEqual(record["protein_chain"], "E")
        self.assertEqual(record["rna_chain"], "A")
        self.assertEqual(record["interface_contacts"], 12)
        self.assertEqual(record["experimental_method"], "x-ray")
        self.assertIsNone(record["backbone_atom_complete"])

    def test_family_annotation_uses_local_assignments(self) -> None:
        record = normalize_record({
            "protein_sequence": "ACDE", "rna_sequence": "ACGU",
            "protein_accession": "P1", "rna_accession": "RNA1.1",
            "source_record_id": "x1", "source_database": "fixture",
            "interaction_evidence": "direct biochemical binding",
        }, "fixture")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            pkey, rkey = record["sequence_key"].split("_r_", 1)
            (root / "protein90_cluster.tsv").write_text(f"{pkey}\t{pkey}\n", encoding="utf-8")
            for field in ("protein40", "protein30"):
                (root / f"{field}_cluster.tsv").write_text(f"{pkey}\t{pkey}\n", encoding="utf-8")
            for field in ("rna90", "rna80"):
                (root / f"{field}_cluster.tsv").write_text(f"r_{rkey}\tr_{rkey}\n", encoding="utf-8")
            (root / "rna_rfam_by_id.tsv").write_text("RNA1\tRF00001\n", encoding="utf-8")
            stats = annotate_family_records([record], root)
            sequence_record = dict(record)
            sequence_record["rna_id"] = "RNA2"
            (root / "rna_rfam_by_seq.tsv").write_text(f"r_{rkey}\tRF00002\n", encoding="utf-8")
            annotate_family_records([sequence_record], root)
        self.assertTrue(record["family_ready"])
        self.assertEqual(record["rfam_family"], "RF00001")
        self.assertEqual(stats["family_ready_rows"], 1)
        self.assertEqual(sequence_record["rfam_family"], "RF00002")


if __name__ == "__main__":
    unittest.main()
