# Final dataset forms

The final train/validation/test tables belong in this directory:

- `train.jsonl` — exactly 100,000 mother samples;
- `validation.jsonl` — experimental-only, with bilateral Protein30/Rfam holdout;
- `test.jsonl` — experimental-only, with bilateral Protein30/Rfam holdout.

These files are generated only after `pri-general validate --strict` passes. They are ignored until a validated release exists.
