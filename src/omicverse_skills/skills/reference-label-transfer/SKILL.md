---
name: omicverse-reference-label-transfer
description: Transfer cell labels from a reference AnnData to a query AnnData with OmicVerse AnnotationRef. Use when converting OmicVerse reference-annotation notebooks into a reusable, triggerable skill, choosing a label-transfer backend such as harmony, scVI, or scanorama, or preparing paired query/reference single-cell data for weighted kNN annotation.
---

# OmicVerse Reference Label Transfer

## Goal

Map labels from a reference `AnnData` onto a query `AnnData` by building a shared integrated space, training a weighted kNN transfer model, and writing predicted labels plus uncertainty scores back to the query object.

Treat `harmony`, `scVI`, and `scanorama` as alternative integration branches. They share the same query/reference contract, so keep them in one skill rather than splitting into thin wrappers.

## Quick Workflow

1. Start from a query `AnnData` and a reference `AnnData` that share gene names.
2. Make sure the reference has a cell-type column such as `celltype` in `adata_ref.obs`.
3. If the data are already log-normalized, recover or restore counts before concatenating.
4. Construct `ov.single.AnnotationRef(adata_query, adata_ref, celltype_key=...)`.
5. Run `preprocess(mode=...)` on the combined object.
6. Choose one integration backend with `train(method=...)`.
7. Transfer labels with `predict(method=...)`.
8. Validate the method-specific prediction keys before plotting or writing outputs.

## Interface Summary

- `AnnotationRef(adata_query, adata_ref, celltype_key='celltype')` checks shared genes, concatenates query and reference, and stores the integrated object on `self.adata_new`.
- `AnnotationRef.preprocess(mode='shiftlog|pearson', n_HVGs=3000, batch_key='integrate_batch')` prepares the concatenated object for transfer.
- `AnnotationRef.train(method='harmony', **kwargs)` computes the integrated latent space used for transfer.
- `AnnotationRef.predict(method='harmony', n_neighbors=15, pred_key=None, uncert_key=None)` transfers labels with weighted kNN.
- The current source also exposes `batch_correction(..., methods='harmony')` with additional methods such as `combat`, `scanorama`, `scVI`, `CellANOVA`, and `Concord`; the notebook uses only the `harmony` and `scVI` transfer path.

Read `references/source-grounding.md` before documenting narrower parameter behavior than the live source supports.

## Branch Selection

- Use `harmony` when you want the notebook's default transfer path and a fast, interpretable integration backend.
- Use `scVI` when you want the scVI latent-space branch and the environment has `scvi-tools`.
- Use `scanorama` when you explicitly want Scanorama integration and its dependencies are available.
- Keep `harmony` as the default unless the user asks for a different integration backend.
- Keep `mode='shiftlog|pearson'` as the default preprocessing choice unless the user wants a different normalization/HVG strategy.

## Input Contract

- Query and reference must share feature names before `AnnotationRef` can concatenate them.
- The reference must carry a cell-type label column in `adata_ref.obs`.
- If the notebook reference was stored with a `feature_name` column, align `adata_ref.var_names` first.
- If `adata_query.X.max()` or `adata_ref.X.max()` suggests the data are already log-normalized, recover counts before transfer.
- `AnnotationRef.preprocess(...)` uses the combined query/reference object, so treat the pair as one integrated preprocessing unit.

## Minimal Execution Patterns

For harmony-based transfer:

```python
import omicverse as ov

ref_anno = ov.single.AnnotationRef(adata_query, adata_ref, celltype_key="celltype")
ref_anno.preprocess(mode="shiftlog|pearson", n_HVGs=3000, batch_key="integrate_batch")
ref_anno.train(method="harmony", n_pcs=50)
ad_pred = ref_anno.predict(method="harmony", n_neighbors=15)
ov.pp.mde(ad_pred, use_rep="X_pca_harmony_anno")
ov.pl.embedding(ad_pred, basis="X_mde", color="harmony_prediction")
```

For scVI-based transfer:

```python
import omicverse as ov

ref_anno = ov.single.AnnotationRef(adata_query, adata_ref, celltype_key="celltype")
ref_anno.preprocess(mode="shiftlog|pearson", n_HVGs=3000, batch_key="integrate_batch")
ref_anno.train(method="scVI", n_layers=2, n_latent=30, gene_likelihood="nb")
ad_pred = ref_anno.predict(method="scVI", n_neighbors=15)
ov.pl.embedding(ad_pred, basis="X_mde", color="scVI_prediction")
```

For prompt-driven label transfer around the reference selection step, keep the notebook's reference selection logic separate from the transfer logic; this skill starts once you already have the paired query/reference objects.

## Validation

After `preprocess(...)`, check:

- `self.adata_new` has `integrate_batch` in `obs`
- `self.adata_new.obsm["X_pca"]` exists
- `self.adata_new.var["highly_variable_features"]` exists

After `train(method="harmony")`, check:

- `self.adata_query.obsm["X_pca_harmony_anno"]` exists
- `self.adata_ref.obsm["X_pca_harmony_anno"]` exists

After `train(method="scVI")`, check:

- `self.adata_query.obsm["X_scVI_anno"]` exists
- `self.adata_ref.obsm["X_scVI_anno"]` exists

After `train(method="scanorama")`, check:

- `self.adata_query.obsm["X_scanorama_anno"]` exists
- `self.adata_ref.obsm["X_scanorama_anno"]` exists

After `predict(...)`, check:

- `harmony_prediction` and `harmony_uncertainty` for the harmony branch
- `scVI_prediction` and `scVI_uncertainty` for the scVI branch
- `scanorama_prediction` and `scanorama_uncertainty` for the Scanorama branch

Before any plotting, check:

- the prediction columns are present in `ad_pred.obs`
- the integrated embedding key you pass to `ov.pp.mde(...)` exists

## Constraints

- Do not bake notebook-specific download URLs, file names, or absolute paths into the reusable workflow.
- Do not assume the reference dataset is already in the correct gene-name convention; align feature names first when needed.
- Do not treat `scVI` as free: it depends on `scvi-tools` and can be materially heavier than the harmony path.
- Do not collapse the notebook into one monolithic "load data and run" step; the reusable spine is query/reference construction, preprocessing, backend selection, prediction, and validation.
- Treat `scanorama` as a source-present branch even though the notebook focuses on harmony and scVI.

## Resource Map

- Read `references/integration-method-selection.md` when choosing between `harmony`, `scVI`, and `scanorama`.
- Read `references/source-notebook-map.md` to trace notebook sections into this reusable skill.
- Read `references/source-grounding.md` for inspected signatures, source-level branch behavior, and compatibility notes.
- Read `references/compatibility.md` when notebook prose, input assumptions, or current runtime behavior diverge.
