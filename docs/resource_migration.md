# Resource migration record

The duplicate `PRI-General_available_resources/` folder and its ZIP were removed after migration. Their former locations contain `.moved.md` pointers.

Migrated local sources:

- `data/cache/sources/npinter_v5/` — official NPInter v5 ZIP, binding-site table, and miRNA table;
- `data/cache/sources/rcsb_pdb/` — PDB candidates and 10,254 experimental-structure assets;
- `data/cache/sources/rnainter/`, `intact/`, `encode_rbns/`, and `rfam/`;
- `data/cache/sources/uniprot/` and `rnacentral/` legacy caches;
- `data/cache/sources/prdbv3/` and `rna3db/` supplemental benchmark/annotation archives.

The old 2.26 GB audit CSV is local-only at `data/cache/sources/legacy_v1/general_stage_audit.csv`; it is not uploaded. Duplicate v1 reports, plots, and derived tables were not copied into v2. The v1 implementation remains recoverable from branch `legacy/v1` and tag `legacy-v1-2026-09-05`.

The official NPInter main download is retained as `interaction_NPInterv5.zip`; its inner table is an evidence/ID source and still needs sequence resolution before it can enter the canonical candidate table.
