# v2 workflow

The pipeline follows one direction:

`raw source -> normalize -> evidence/construct QC -> experimental interface QC -> family assignment -> deterministic selection -> isolated split -> external Boltz result validation -> release`

The code intentionally stops at a draft release when an upstream contract is missing. It does not convert a database row count into a biological sample count, and it does not treat a prediction as experimental evidence.

## Current source disposition

- `npinter_bindingsite`, `npinter_mirna`, and the RCSB candidate table are the current sequence-ready candidate inputs.
- `npinter_main` is a valid local ZIP evidence/ID source, but its identifiers still require sequence resolution; it is not imported as guessed rows.
- RNAcentral current ID mapping and Rfam annotations were downloaded, gzip-verified, and used for local family annotation. Files over 1 GB remain local-only.
- RNApedia is a verified 35.5 GB Zenodo benchmark archive and remains local-only; it is not mixed into training candidates.
- RNAInter, IntAct, and ENCODE/RBNS are cataloged in the same manifest. Their current local materials do not provide a complete sequence-resolved canonical table, so they remain provenance/annotation inputs until that requirement is met.
- POSTAR3 has a public web endpoint but no verified machine-readable dump was found in this pass.
- Boltz-2/MSA are deliberately external to this project.
