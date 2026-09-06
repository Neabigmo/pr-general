from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path
from typing import Iterable

from Bio.PDB import FastMMCIFParser, NeighborSearch
from Bio.PDB.Polypeptide import is_aa, is_nucleic


def _atom_name(atom) -> str:
    return atom.get_name().strip().replace("*", "'")


def _chain_atoms(chain, predicate):
    residues = {}
    for residue in chain:
        if not predicate(residue):
            continue
        atoms = [atom for atom in residue if not _atom_name(atom).startswith("H")]
        if atoms:
            residues[residue.get_id()] = (residue, atoms)
    return residues


def _pair_metrics(model, protein_chain_id: str, rna_chain_id: str, cutoff: float) -> dict:
    chains = {chain.id: chain for chain in model}
    protein_chain = chains.get(protein_chain_id)
    rna_chain = chains.get(rna_chain_id)
    if protein_chain is None or rna_chain is None:
        return {"backbone_atom_complete": False, "interface_qc_reason": "missing_chain"}
    protein = _chain_atoms(protein_chain, lambda residue: is_aa(residue, standard=False))
    rna = _chain_atoms(rna_chain, lambda residue: is_nucleic(residue, standard=False))
    if not protein or not rna:
        return {"backbone_atom_complete": False, "interface_qc_reason": "missing_protein_or_rna_atoms"}
    rna_atoms = [atom for _, atoms in rna.values() for atom in atoms]
    nearby = NeighborSearch(rna_atoms)
    protein_contacts = set()
    rna_contacts = set()
    contact_pairs = 0
    for protein_id, (_, atoms) in protein.items():
        for atom in atoms:
            for rna_atom in nearby.search(atom.coord, cutoff, level="A"):
                contact_pairs += 1
                protein_contacts.add(protein_id)
                rna_contacts.add(rna_atom.get_parent().get_id())
    protein_complete = all({"N", "CA", "C", "O"}.issubset({_atom_name(atom) for atom in protein[residue_id][1]}) for residue_id in protein_contacts)
    rna_complete = all(
        {"C3'", "C4'", "C5'"}.issubset({_atom_name(atom) for atom in rna[residue_id][1]})
        and bool({"P", "O3'", "O5'"} & {_atom_name(atom) for atom in rna[residue_id][1]})
        for residue_id in rna_contacts
    )
    complete = protein_complete and rna_complete
    return {
        "interface_contacts": contact_pairs,
        "contact_protein_residues": len(protein_contacts),
        "contact_rna_nucleotides": len(rna_contacts),
        "backbone_atom_complete": complete,
        "interface_qc_reason": "qualified_metrics" if complete else "incomplete_backbone_atoms",
    }


def scan_partition(rows: Iterable[dict], structure_root: str | Path, output_path: str | Path, partition: int, partitions: int, cutoff: float = 5.0) -> dict[str, int]:
    """Scan a disjoint PDB partition; each worker writes its own JSONL output."""
    grouped: dict[str, dict[tuple[str, str], list[str]]] = defaultdict(lambda: defaultdict(list))
    for row in rows:
        pdb_id = str(row.get("structure_id") or "").upper()
        if not pdb_id or sum(ord(char) for char in pdb_id) % partitions != partition:
            continue
        protein_chain = str(row.get("protein_chain") or "")
        rna_chain = str(row.get("rna_chain") or "")
        if protein_chain and rna_chain:
            grouped[pdb_id][(protein_chain, rna_chain)].append(str(row["sample_id"]))
    parser = FastMMCIFParser(QUIET=True)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    scanned = failed = written = 0
    with output.open("w", encoding="utf-8") as handle:
        for pdb_id, pairs in sorted(grouped.items()):
            cif = Path(structure_root) / f"{pdb_id}-assembly1.cif"
            try:
                structure = parser.get_structure(pdb_id, cif)
                model = next(structure.get_models())
                for (protein_chain, rna_chain), sample_ids in pairs.items():
                    metrics = _pair_metrics(model, protein_chain, rna_chain, cutoff)
                    for sample_id in sample_ids:
                        handle.write(json.dumps({"sample_id": sample_id, **metrics}, sort_keys=True) + "\n")
                        written += 1
                scanned += 1
                if scanned % 25 == 0:
                    handle.flush()
            except Exception as exc:  # retain an auditable failure instead of guessing metrics
                failed += 1
                for sample_ids in pairs.values():
                    for sample_id in sample_ids:
                        handle.write(json.dumps({"sample_id": sample_id, "backbone_atom_complete": False, "interface_qc_reason": f"parse_error:{type(exc).__name__}"}, sort_keys=True) + "\n")
                        written += 1
    return {"pdb_groups": len(grouped), "pdb_scanned": scanned, "pdb_failed": failed, "rows_written": written}


def apply_structure_qc(rows: list[dict], qc_path: str | Path) -> int:
    metrics = {}
    path = Path(qc_path)
    if not path.exists():
        return 0
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                item = json.loads(line)
                metrics[item["sample_id"]] = item
    applied = 0
    for row in rows:
        item = metrics.get(row.get("sample_id"))
        if item is None:
            continue
        row.update({key: value for key, value in item.items() if key != "sample_id"})
        applied += 1
    return applied
