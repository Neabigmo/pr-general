from __future__ import annotations

from typing import Any

from .constants import clean_sequence, nonempty, pair_digest, sequence_digest
from .evidence import infer_evidence_type


def _first(raw: dict, *names: str) -> Any:
    for name in names:
        value = raw.get(name)
        if value is not None and str(value).strip() != "":
            return value
    return None


def normalize_record(raw: dict, source_id: str) -> dict:
    protein = clean_sequence(_first(raw, "protein_sequence", "protein_seq", "protein"), rna=False)
    rna = clean_sequence(_first(raw, "rna_sequence", "rna_seq", "rna"), rna=True)
    protein_id = nonempty(_first(raw, "protein_accession", "protein_id", "uniprot_id", "tarID"))
    rna_id = nonempty(_first(raw, "rna_accession", "rna_id", "rnacentral_id", "ncID"))
    source_record_id = nonempty(_first(raw, "source_record_id", "record_id", "sample_id", "interID"))
    structure_id = nonempty(_first(raw, "pdb_id", "structure_id"))
    start = _first(raw, "binding_site_start", "window_start", "site_start")
    end = _first(raw, "binding_site_end", "window_end", "site_end")
    evidence_type = infer_evidence_type(raw)
    pair_identity = "|".join((protein_id or protein, rna_id or rna))
    sample_id = f"{source_id}:{source_record_id}" if source_record_id else f"{source_id}:{pair_digest(protein, rna)[:16]}"
    return {
        "sample_id": sample_id,
        "biological_pair_id": f"bp_{pair_digest(pair_identity, '')[:24]}",
        "protein_id": protein_id,
        "rna_id": rna_id,
        "protein_sequence": protein,
        "rna_sequence": rna,
        "protein_length": len(protein),
        "rna_length": len(rna),
        "protein_gene_id": nonempty(_first(raw, "protein_gene_id", "gene_id")),
        "rna_gene_id": nonempty(_first(raw, "rna_gene_id", "gene_id_rna")),
        "organism": nonempty(raw.get("organism")),
        "rna_class": nonempty(_first(raw, "rna_class", "system_subtype")),
        "evidence_type": evidence_type,
        "evidence_id": nonempty(_first(raw, "evidence_id", "source_record_id", "record_id")),
        "evidence_support": _first(raw, "evidence_support", "confidence", "score"),
        "source_database": nonempty(_first(raw, "source_database", "source", "database", "datasource")) or source_id,
        "source_record_id": source_record_id,
        "structure_id": structure_id,
        "protein_chain": nonempty(_first(raw, "protein_chain", "chain_protein")),
        "rna_chain": nonempty(_first(raw, "rna_chain", "chain_rna")),
        "site_mode": "site_window" if start is not None and end is not None else "full_construct",
        "binding_site_start": start,
        "binding_site_end": end,
        "experimental_method": nonempty(_first(raw, "experimental_method", "method")),
        "resolution_angstrom": _first(raw, "resolution_angstrom", "resolution"),
        "interface_contacts": _first(raw, "interface_contacts", "n_contacts", "contact_pairs"),
        "contact_protein_residues": _first(raw, "contact_protein_residues", "n_protein_residues"),
        "contact_rna_nucleotides": _first(raw, "contact_rna_nucleotides", "n_rna_nucleotides"),
        "backbone_atom_complete": _first(raw, "backbone_atom_complete", "backbone_atoms_complete"),
        "quality_tier": None,
        "protein90": nonempty(raw.get("protein90")),
        "protein40": nonempty(raw.get("protein40")),
        "protein30": nonempty(raw.get("protein30")),
        "rna90": nonempty(raw.get("rna90")),
        "rna80": nonempty(raw.get("rna80")),
        "rfam_family": nonempty(_first(raw, "rfam_family", "rfam")),
        "rfam_clan": nonempty(raw.get("rfam_clan")),
        "split": None,
        "sequence_key": f"p_{sequence_digest(protein)[:24]}_r_{sequence_digest(rna)[:24]}",
        "source_id": source_id,
    }
