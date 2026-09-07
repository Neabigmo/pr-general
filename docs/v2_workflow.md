# v2 workflow

The pipeline follows one direction:

`raw source -> normalize -> evidence/construct QC -> experimental interface QC -> family assignment -> deterministic selection -> isolated split -> external Boltz result validation -> release`

The code intentionally stops at a draft release when an upstream contract is missing. It does not convert a database row count into a biological sample count, and it does not treat a prediction as experimental evidence.

## Current source disposition

- `npinter_bindingsite`, `npinter_mirna`, the resolved NPInter main subset, and the RCSB candidate table are the current sequence-ready candidate inputs.
- Raw source caches retain the source-level provenance, but only one evidence-prioritized representative per exact sequence pair enters the canonical mother-sample pool.
- `npinter_main` keeps the raw ZIP as its source-of-record and uses a local resolved subset built only from exact miRBase hairpins plus exact UniProt accessions or unique RefSeq CDS translations for exact NM accessions; unresolved identifiers remain provenance-only.
- RNAcentral current ID mapping and Rfam annotations were downloaded, gzip-verified, and used for local family annotation. Files over 1 GB remain local-only.
- Training readiness requires Protein30 and RNA80. Rfam is preferred for RNA family isolation; when it is absent, the split logic uses RNA80 and does not consume an Rfam-family cap or invent an annotation.
- The current local snapshot produces 33,666 cap-valid, exact-sequence mother samples after evidence-prioritized collapse and strict holdout reservation, before external structure prediction; the configured 100,000 target remains blocked by the exact-protein cap and available sequence diversity.
- RNApedia is a verified 35.5 GB Zenodo benchmark archive and remains local-only; it is not mixed into training candidates.
- RNAInter, IntAct, and ENCODE/RBNS are cataloged in the same manifest. Their current local materials do not provide a complete sequence-resolved canonical table, so they remain provenance/annotation inputs until that requirement is met.
- POSTAR3's web query endpoint responds, but no verified machine-readable bulk dump was found in this pass; ad hoc query results are not imported as a bulk source.
- Boltz-2/MSA are deliberately external to this project.
