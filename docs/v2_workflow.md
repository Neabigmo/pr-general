# v2 workflow

The pipeline follows one direction:

`raw source -> normalize -> evidence/construct QC -> experimental interface QC -> family assignment -> deterministic selection -> isolated split -> external Boltz result validation -> release`

The code intentionally stops at a draft release when an upstream contract is missing. It does not convert a database row count into a biological sample count, and it does not treat a prediction as experimental evidence.

## Current resource gaps

- `npinter_main` is present locally but corrupt: it has a Parquet header without a footer.
- RNAcentral current ID mapping and Rfam annotations are available from the official FTP but are large downloads and are recorded as pending.
- RNApedia is available from its official download page, but its complete archive is about 35.5 GB and overlaps the local PDB structural cache.
- POSTAR3 has a public web endpoint but no verified machine-readable dump was found in this pass.
- Boltz-2/MSA are deliberately external to this project.
