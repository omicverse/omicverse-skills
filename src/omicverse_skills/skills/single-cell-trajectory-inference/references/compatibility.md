# Compatibility

## Upstream Reuse

This notebook starts with a preprocessing block that is already a standalone skill boundary. Do not duplicate that execution spine here; use the preprocessing skill first when the object is not already trajectory-ready.

## Dependency-Sensitive Branches

### Slingshot

Current OmicVerse behavior:

- `TrajInfer.inference(method='slingshot', ...)` imports `Slingshot` from `omicverse/single/_pyslingshot.py`
- `Slingshot.fit(...)` imports `pcurvepy2.PrincipalCurve`

Consequence:

- the branch is unavailable unless `pcurvepy2` is installed

### Palantir Gene Trends

Current OmicVerse behavior:

- `palantir_cal_gene_trends(...)` delegates to `compute_gene_trends(...)`
- that helper imports `mellon`

Consequence:

- the trend branch is unavailable unless `mellon` is installed

## Palantir Wrapper Limitation

The wrapper in `omicverse/single/_traj.py` calls `run_magic_imputation(self.adata)` directly inside the Palantir branch and does not expose the helper's `n_jobs` parameter on the public `TrajInfer.inference(...)` surface.

Consequence:

- constrained local smoke harnesses may need a single-process patch when validating the branch under restrictive process limits
- keep that workaround in the validation harness only; do not promote it into the reusable workflow

## Notebook-Specific Labels

The notebook uses pancreas-specific states such as `Ductal`, `Alpha`, `Beta`, `Delta`, and `Epsilon`.

Treat them as examples only. Another compatible dataset should supply its own biological labels for:

- origin
- optional terminal states
- marker genes used in follow-up plots
