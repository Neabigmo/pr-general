# PRI-General v2

PRI-General v2 is a clean rebuild of the protein–RNA dataset pipeline. The old implementation is frozen on branch `legacy/v1` and tag `legacy-v1-2026-09-05`; this branch contains only the v2 workflow.

## Scope

- one biological protein–RNA pair per mother sample;
- experimental evidence is separated from computational evidence;
- structure QC is pair-level and requires explicit interface metrics;
- Protein90/40/30 and RNA90/80/Rfam assignments are inputs to selection and split isolation;
- validation and test are experimental-only and isolated by biological pair, Protein30, and Rfam;
- Boltz-2 execution and MSA preparation are external to this project, as requested.

## Run

```bash
python -m pip install -e .
python -m pri_general build --config config/pri_general_v2.yaml --release v2.0.0
python -m pri_general validate --config config/pri_general_v2.yaml
python -m pri_general report --config config/pri_general_v2.yaml --release v2.0.0
```

`build` writes local intermediate data to `data/cache/build/` and never places raw data or generated tables in Git. A build without family assignments and external Boltz results remains `DRAFT`; it must not be presented as the final dataset.

## Important locations

- source inventory: `data/manifests/sources.yaml`;
- local raw/cache data: `data/cache/sources/` and `data/cache/build/`;
- final dataset tables/forms: `releases/v2.0.0/forms/`;
- prediction submission template and completed submission: `releases/v2.0.0/submission/`;
- release metadata and audit reports: `releases/v2.0.0/`.

Files over 1 GB are local-only and are represented by small explanation files when a repository location is needed.
