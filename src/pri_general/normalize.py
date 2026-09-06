from __future__ import annotations

import math
import re
from typing import Any

from .constants import clean_sequence, nonempty, pair_digest, sequence_digest
from .evidence import infer_evidence_type


def _first(raw: dict, *names: str) -> Any:
    for name in names:
        value = raw.get(name)
        if value is None:
            continue
        try:
            if math.isnan(value):
                continue
        except TypeError:
            pass
        if str(value).strip() != "":
            return value
    return None


def _note_value(notes: str, name: str) -> str | None:
    match = re.search(rf"(?:^|;)\s*{re.escape(name)}\s*=\s*([^;]+)", notes)
    return match.group(1).strip() if match else None


def _structure_fields(raw: dict) -> dict[str, Any]:
    notes = str(raw.get("notes") or "")
    evidence = str(raw.get("interaction_evidence") or "")
    protein_accession = str(_first(raw, "protein_accession") or "")
    rna_accession = str(_first(raw, "rna_accession") or "")
    protein_chain = _first(raw, "protein_chain", "chain_protein") or _note_value(notes, "protein_chain")
    rna_chain = _first(raw, "rna_chain", "chain_rna") or _note_value(notes, "rna_chain")
    if not protein_chain and ":" in protein_accession:
        protein_chain = protein_accession.rsplit(":", 1)[1]
    if not rna_chain and ":" in rna_accession:
        rna_chain = rna_accession.rsplit(":", 1)[1]
    resolution = _first(raw, "resolution_angstrom", "resolution")
    if resolution is None:
        match = re.search(r"resolution\s*=\s*\[?([0-9]+(?:\.[0-9]+)?)", notes, re.I)
        resolution = float(match.group(1)) if match else None
    contacts = _first(raw, "interface_contacts", "n_contacts", "contact_pairs")
    if contacts is None:
        match = re.search(r"(\d+)\s+contact atom pairs", evidence, re.I)
        contacts = int(match.group(1)) if match else None
    method = _first(raw, "experimental_method", "method")
    if method is None:
        method = _note_value(notes, "method")
    method = {"X": "x-ray", "N": "nmr", "E": "electron microscopy", "C": "cryo-em"}.get(
        str(method).strip().upper(), method
    )
    return {
        "protein_chain": protein_chain,
        "rna_chain": rna_chain,
        "resolution_angstrom": resolution,
        "interface_contacts": contacts,
        "experimental_method": method,
    }


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
    sample_suffix = pair_digest(protein, rna)[:16]
    site_suffix = f":site-{start}-{end}" if start is not None and end is not None else ""
    sample_id = (
        f"{source_id}:{source_record_id}:{sample_suffix}{site_suffix}"
        if source_record_id
        else f"{source_id}:{sample_suffix}{site_suffix}"
    )
    structure = _structure_fields(raw)
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
        "protein_chain": nonempty(structure["protein_chain"]),
        "rna_chain": nonempty(structure["rna_chain"]),
        "site_mode": "site_window" if start is not None and end is not None else "full_construct",
        "binding_site_start": start,
        "binding_site_end": end,
        "experimental_method": nonempty(structure["experimental_method"]),
        "resolution_angstrom": structure["resolution_angstrom"],
        "interface_contacts": structure["interface_contacts"],
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
