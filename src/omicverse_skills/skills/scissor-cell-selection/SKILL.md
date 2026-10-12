---
name: scissor-cell-selection
description: Select bulk-phenotype-associated single cells with the independent scissor_py kernel or ov.single.Scissor, using explicit expression layers, cell graphs and phenotype coding. Use for Gaussian, binary or right-censored survival Scissor analyses and their reproducibility checks.
---

# Scissor cell selection

Use the user's supplied data and analysis choices. The backend import is
`scissor_py`; the unrelated `pyscissor` distribution is not a substitute. The
OmicVerse integration is optional. Check the installed version and whether
`ov.single.Scissor` exists before using it; the adapter PR must be merged or its reviewed revision installed before use. Do not install packages, download datasets or publish results
as an implicit side effect of a library call.

## Backend and version boundary

The distribution is `scissor-py` (Python >=3.10). When installation is requested,
use this verified organization commit:

```sh
python -m pip install "git+https://github.com/omicverse/py-Scissors.git@cc630ae828149202b553fb7dfdaeda2d0b94d1d7"
```

Check the active interpreter and record the installed backend and OmicVerse
versions. An earlier checkout may also report `0.1.0.dev5`; when replacement
is requested, add `--force-reinstall --no-deps` to the pinned backend command
to replace that revision without reinstalling the scientific dependencies. Installing the backend does not add `ov.single.Scissor` to an older
OmicVerse version. For a missing adapter, use the independent backend or the
[reviewed OmicVerse contribution](https://github.com/omicverse/omicverse/pull/988) explicitly.
The validated adapter revision is `b04249097f35ab261fcbc6e30becfb1d85cb046b`
in `hurry060215-tech/omicverse`; it is not yet a merged upstream release.

## Scientific input decisions

- Bulk is a pandas DataFrame, genes × samples. Single-cell AnnData is cells ×
  genes. Require unique, nonmissing IDs; phenotype and supplied fold IDs must
  match the bulk samples exactly. Do not intersect away unmatched samples or
  silently resolve duplicate genes.
- Require the intended comparable expression layer and an existing cell graph.
  In `ov.single.Scissor`, `layer=None` explicitly selects X. Missing graphs need
  a chosen upstream preprocessing/neighbors workflow, not a silently created
  replacement. Count/log normalization, HVG selection and batch correction are
  outside this kernel. Different graphs, including Scanpy versus Seurat SNN,
  can change selected cells.
- Gaussian uses an indexed numeric Series. Binomial uses an indexed two-class
  Series; require `positive_label` for strings, and retain the positive/negative
  mapping. Numeric 0/1 defaults to positive 1. Cox uses exactly time/status
  columns; 1 means event, 0 means right-censoring. Time ties use Breslow.
- Explain numerical semantics before choosing them. `official` preserves
  defined behavior of fixed Scissor commit 311560a1 and rejects unsafe
  Gaussian/Cox graph indexing. `corrected` explicitly changes graph alignment,
  Gaussian curvature, stopping rules and Cox lambda bookkeeping. Never silently
  switch semantics after an error. Read the installed project's
  [semantics documentation](https://github.com/omicverse/py-Scissors/blob/cc630ae828149202b553fb7dfdaeda2d0b94d1d7/docs/SEMANTICS.md) for the compatibility boundary.

## Fit and retain evidence

```python
model = ov.single.Scissor(
    adata, bulk=bulk, phenotype=phenotype, family="binomial",
    layer="lognorm", graph_key="connectivities", positive_label="mutant",
    semantics="corrected", key_added="scissor", overwrite=False,
)
result = model.run(alpha=None, cutoff=.2, random_state=123, foldid=foldid)
```

The first successful alpha with selected fraction strictly below cutoff wins.
Each alpha uses CV lambda.min followed by a cold single-lambda refit. Equal R
and Python seeds do not guarantee equal folds: preserve actual foldid, graph,
input ordering, alpha and lambda path for comparisons. Gaussian CV uses held-out
MSE, binary CV deviance, and Cox event-weighted full-minus-training partial
likelihood, not C-index. Nonconverged candidates cannot win; all failures must
raise. If cutoff is not reached, report that status rather than implying sparse
selection succeeded.

The optional Numba engine is explicit and currently only supports corrected
Gaussian. Default NumPy requires no Numba or R. Do not change solver settings,
fold count, lambda count or sampled data merely to make a validation pass.
Tighter convergence is a separate recorded experiment, not permission to relax
acceptance tolerances or hide default-mode discrepancies.

## Outputs and interpretation

Successful annotation writes obs scissor/scissor_coef and uns scissor, or the
chosen key namespace. Existing keys require explicit overwrite; views/backed
objects need an explicit copy/to_memory before writing. Core fitting can remain
read-only through `scissor_py.scissor_from_anndata(..., writeback=False)`.

AnnData displays Scissor+/Scissor-/background. Returned core results retain
Background for the zero coefficient label. Selection uses strict signs; near-zero
coefficients are diagnostic, not silently relabelled. Scissor+ points toward
the positive binary class, greater continuous outcome or higher Cox log hazard.
These are model associations, not causal or clinical conclusions.

Inspect cutoff_reached, failed candidates, convergence, near-zero cells and
metadata. Save the standard-JSON audit and, if requested, use the separate
`write_h5ad_atomic` helper for no-half-file output. The audit excludes dense
matrices and full-path coefficients; preserve those separately when a full
rerun or model reload is needed. H5AD JSON fields require json.loads and retain
ID/class types and explicit nonfinite-diagnostic tags. See the [AnnData serialization contract](https://github.com/omicverse/py-Scissors/blob/cc630ae828149202b553fb7dfdaeda2d0b94d1d7/docs/ANNDATA.md).
