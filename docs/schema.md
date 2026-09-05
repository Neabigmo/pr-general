# Canonical record

The canonical record is defined in `src/pri_general/constants.py`. Identity is represented at three levels:

- `sample_id`: source-row-level reproducible identifier;
- `biological_pair_id`: one protein/RNA biological pair for grouping and split isolation;
- `sequence_key`: exact normalized sequence identity for deduplication and caps.

The content keys are necessary because source row IDs and database releases are not stable biological identities. They are not release gates or a replacement for ordinary tests.
