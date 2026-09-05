from __future__ import annotations

import hashlib
import re

PROTEIN_ALPHABET = set("ACDEFGHIKLMNPQRSTVWYBXZJUO")
RNA_ALPHABET = set("ACGUNRYSWKMBDHV")

EVIDENCE_TYPES = (
    "experimental_structure",
    "direct_biochemical",
    "site_resolved_binding",
    "physical_association",
    "computational_prediction",
)

CANONICAL_FIELDS = (
    "sample_id", "biological_pair_id", "protein_id", "rna_id",
    "protein_sequence", "rna_sequence", "protein_length", "rna_length",
    "protein_gene_id", "rna_gene_id", "organism", "rna_class",
    "evidence_type", "evidence_id", "evidence_support", "source_database",
    "source_record_id", "structure_id", "protein_chain", "rna_chain",
    "site_mode", "binding_site_start", "binding_site_end",
    "experimental_method", "resolution_angstrom", "interface_contacts",
    "contact_protein_residues", "contact_rna_nucleotides",
    "backbone_atom_complete", "quality_tier", "protein90", "protein40",
    "protein30", "rna90", "rna80", "rfam_family", "rfam_clan", "split",
)


def clean_sequence(value: object, *, rna: bool) -> str:
    if value is None:
        return ""
    text = re.sub(r"\s+", "", str(value)).upper()
    return text.replace("T", "U") if rna else text


def sequence_digest(sequence: str) -> str:
    # ponytail: a stable content key is necessary because database row IDs and
    # release/version labels are not stable biological identities.
    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()


def pair_digest(left: str, right: str) -> str:
    return sequence_digest(f"{left}\x1f{right}")


def nonempty(value: object) -> str | None:
    text = "" if value is None else str(value).strip()
    return text or None
