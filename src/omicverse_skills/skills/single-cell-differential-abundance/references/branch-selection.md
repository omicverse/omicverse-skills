# Branch Selection

## Backend Choice

- Choose `sccoda` for Bayesian compositional inference with posterior sampling and optional FDR tuning after fitting.
- Choose `milopy` for the native Python Milo-like neighborhood path.
- Choose `milo` only when the pertpy Milo branch is explicitly required.

## Required Inputs Per Branch

- `sccoda`: condition labels, cell-type labels, and a sample identifier.
- `milopy`: condition labels, cell-type labels, sample identifier, and an embedding already stored in `obsm`.
- `milo`: the same practical requirements as `milopy`, even though constructor validation is looser.

## Result Interpretation

- Use `mix_threshold` only for the Milo-family result tables.
- Apply posterior-specific follow-up such as FDR adjustment only after the `sccoda` model has completed.

## Plotting

- Add histogram, MA-style, or beeswarm plots only after the abundance result table exists.
- Treat plot aesthetics as request-specific reporting rather than part of the reusable skill contract.
