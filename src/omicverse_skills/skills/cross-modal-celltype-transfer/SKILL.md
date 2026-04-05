---
name: omicverse-cross-modal-celltype-transfer
description: Transfer cell-type labels from a reference AnnData to a query AnnData with OmicVerse weighted KNN over a shared embedding. Use when converting OmicVerse cross-modal annotation notebooks into a reusable skill, when labeling an ATAC query from an RNA reference, or when you need the weighted_knn_trainer / weighted_knn_transfer workflow plus optional visualization.
---

# OmicVerse Cross-Modal Cell Type Transfer

## Goal

Transfer cell labels from a reference `AnnData` to a query `AnnData` using a shared latent space such as `X_glue`, then write predicted labels and uncertainty scores back to the query object. Keep this as one skill: the notebook's plots are optional presentation on top of the same transfer job, not a separate capability.

## Quick Workflow

1. Confirm the query and reference expose the same embedding basis in `.obsm`, or fall back to `X`.
2. Confirm the reference has the target label column in `.obs`.
3. Train `ov.utils.weighted_knn_trainer(...)` on the reference embedding.
4. Transfer labels with `ov.utils.weighted_knn_transfer(...)`.
5. Store the first label column and uncertainty column in `query_adata.obs`.
6. Plot on an existing basis for smoke validation; use `mde` only if `pymde` is available.

## Interface Summary

- `weighted_knn_trainer(train_adata, train_adata_emb, n_neighbors=50)` builds the reference KNN model.
- `weighted_knn_transfer(query_adata, query_adata_emb, ref_adata_obs, label_keys, knn_model, threshold=1, pred_unknown=False, mode='package')` returns prediction and uncertainty DataFrames.
- `embedding(adata, basis, ..., show=None, save=None, return_fig=None, ...)` plots any existing basis.
- `mde(data, device=None, **kwargs)` is optional; it depends on `pymde` and is only for visualization.

## Branch Selection

- Use `mode='package'` for transfer. The current source only implements this branch; `mode='paper'` is documented in the docstring but raises in the live code.
- Leave `pred_unknown=False` for the default label-assignment path.
- Set `pred_unknown=True` only when you want threshold-based `Unknown` labels.
- Use `mde` only as an optional visual branch. If `pymde` is missing, skip it and plot on `X_glue` or another precomputed basis instead.

## Input Contract

- `train_adata` and `query_adata` must expose the same latent basis name in `.obsm`, or use `X`.
- `ref_adata_obs[label_keys]` must exist; the function selects columns whose names start with `label_keys`.
- If you want one label column, make the prefix unique.
- The notebook's `transf_celltype` and `transf_celltype_unc` are the reuse convention for storing outputs on the query object.

## Minimal Execution Patterns

Core transfer:

```python
knn = ov.utils.weighted_knn_trainer(train_adata=rna, train_adata_emb="X_glue", n_neighbors=15)
labels, uncert = ov.utils.weighted_knn_transfer(
    query_adata=atac,
    query_adata_emb="X_glue",
    ref_adata_obs=rna.obs,
    label_keys="major_celltype",
    knn_model=knn,
    mode="package",
)
atac.obs["transf_celltype"] = labels.loc[atac.obs.index, "major_celltype"]
atac.obs["transf_celltype_unc"] = uncert.loc[atac.obs.index, "major_celltype"]
```

Optional smoke plot:

```python
ov.utils.embedding(
    atac,
    basis="X_glue",
    color="transf_celltype",
    show=False,
)
```

## Validation

- `labels` and `uncert` must both be DataFrames indexed by `query_adata.obs_names`.
- The expected label column should appear in both outputs.
- `query_adata.obs["transf_celltype"]` and `query_adata.obs["transf_celltype_unc"]` should exist after writeback.
- If you call `mde`, verify `pymde` is installed first; otherwise skip it.
- For plots, confirm the basis you pass actually exists in `.obsm`.

## Resource Map

- Read `references/integration-method-selection.md` for skill boundary and branch selection.
- Read `references/source-notebook-map.md` to map notebook cells to the reusable workflow.
- Read `references/source-grounding.md` for inspected signatures and source behavior.
- Read `references/compatibility.md` for known runtime caveats.
