---
name: omicverse-single-cell-batch-integration
description: Run OmicVerse single-cell batch integration as a reusable, triggerable skill after preprocessing is already complete. Use when choosing a batch correction backend such as harmony, combat, scanorama, scVI, CellANOVA, or Concord on a preprocessed AnnData with batch labels, or when benchmarking integrated embeddings from one of those backends.
---

# OmicVerse Single-Cell Batch Integration

## Goal

Turn the notebook's post-preprocessing batch integration workflow into one reusable skill. Assume QC, normalization, HVG selection, and base PCA preparation are already handled upstream; this skill starts from a preprocessed `AnnData` with a batch column and focuses on integration backend selection, optional embedding visualization, and optional integration benchmarking.

## Quick Workflow

1. Confirm that preprocessing is already done and that `adata.obs` contains the batch labels.
2. Choose an integration backend explicitly with `methods=...`.
3. Validate the backend-specific prerequisites before calling `ov.single.batch_correction(...)`.
4. Check which integrated embedding key was written to `adata.obsm`, and keep model/object returns when the backend returns more than `adata`.
5. Run optional embedding visualization only after the integrated embedding exists.
6. Run benchmarking only when at least one non-integrated reference embedding and one or more integrated embeddings are present.

## Interface Summary

- `ov.single.batch_correction(adata, batch_key, use_rep='scaled|original|X_pca', methods='harmony', n_pcs=50, **kwargs)` dispatches across the notebook's integration backends.
- `methods='harmony'` writes both `X_pca_harmony` and `X_harmony`.
- `methods='combat'` runs ComBat on a copy, rescales, recomputes PCA, and writes `X_combat`.
- `methods='scanorama'` splits by categorical batch, runs Scanorama per batch list, concatenates the integrated coordinates, and writes `X_scanorama`.
- `methods='scVI'` calls `scvi.model.SCVI.setup_anndata(..., layer='counts', batch_key=batch_key)`, trains a model, writes `X_scVI`, and returns the trained model object instead of returning `adata`.
- `methods='CellANOVA'` expects a `control_dict`, computes `denoised`, runs PCA on that layer, and writes `X_cellanova`.
- `methods='Concord'` or `methods='concord'` constructs a `concord.Concord` object, runs `fit_transform(output_key='X_concord')`, writes `X_concord`, and returns the Concord object.
- `ov.pp.mde(...)` can be used for optional low-dimensional visualization once an embedding key exists.
- `scib_metrics.benchmark.Benchmarker(...)` benchmarks one or more embedding keys against batch and label annotations.

## Stage Selection

- Keep this as one skill because the notebook's integration and benchmarking stages share one tight input contract: a preprocessed `AnnData` with batch labels and integrated embedding outputs.
- Use the integration stage when the user wants a corrected latent space or wants to compare backends.
- Use the visualization stage only after a concrete integrated embedding key exists.
- Use the benchmarking stage only after integration outputs exist; it is downstream evaluation, not an independent upstream job.
- Do not re-run generic preprocessing here; reuse the existing single-cell preprocessing skill first when needed.

## Backend Selection

- Use `harmony` as the default CPU-friendly backend when you already have PCA and want the notebook's primary integration path.
- Use `combat` when you want a simple expression-level correction path and can accept a recomputed PCA output.
- Use `scanorama` when pairwise manifold stitching is explicitly requested and its extra dependencies are available.
- Use `scVI` when the environment has `scvi-tools`, raw counts are present in `adata.layers['counts']`, and a heavier latent-variable model is acceptable.
- Use `CellANOVA` when you have a defensible control pool and explicitly want signal recovery with `control_dict`.
- Use `Concord` when the environment has Concord installed and the user wants its self-supervised latent space.

## Input Contract

- Start from a preprocessed `AnnData` that already has a batch column in `adata.obs`.
- Preserve raw counts in `adata.layers['counts']` before `scVI`.
- Keep or create a PCA representation before `harmony` when you want to control `use_rep` explicitly.
- Convert the batch column to categorical before `scanorama` so the per-batch split is stable.
- Ensure `control_dict` is provided for `CellANOVA` and that it names valid batch groups.
- Keep a non-integrated reference embedding such as `X_pca` if benchmarking will follow.

## Minimal Execution Patterns

```python
import omicverse as ov

ov.single.batch_correction(
    adata,
    batch_key="batch",
    methods="harmony",
    n_pcs=50,
    use_gpu=False,
)
```

```python
ov.single.batch_correction(
    adata,
    batch_key="batch",
    methods="combat",
    n_pcs=50,
)
```

```python
model = ov.single.batch_correction(
    adata,
    batch_key="batch",
    methods="scVI",
    n_layers=2,
    n_latent=30,
    gene_likelihood="nb",
)
```

```python
control_dict = {
    "pool1": ["batch_a", "batch_b"],
}

ov.single.batch_correction(
    adata,
    batch_key="batch",
    methods="CellANOVA",
    n_pcs=50,
    control_dict=control_dict,
)
```

```python
from scib_metrics.benchmark import Benchmarker

bm = Benchmarker(
    adata,
    batch_key="batch",
    label_key="cell_type",
    embedding_obsm_keys=["X_pca", "X_pca_harmony", "X_combat"],
    pre_integrated_embedding_obsm_key="X_pca",
    n_jobs=1,
    progress_bar=False,
)
bm.benchmark()
results = bm.get_results(min_max_scale=False)
```

## Constraints

- Do not treat this as a preprocessing skill.
- Do not leave `methods` implicit.
- Do not assume every backend returns `adata`; keep the returned object for `scVI` and `Concord`.
- Do not claim that all backends are equally cheap; `scVI` and `Concord` are materially heavier than `harmony` or `combat`.
- Do not benchmark embeddings that were never written to `adata.obsm`.
- Keep smoke commands shell-agnostic.

## Validation

- After `harmony`, check `X_pca_harmony` and `X_harmony`.
- After `combat`, check `X_combat`.
- After `scanorama`, check `X_scanorama`.
- After `scVI`, check that a model object was returned and that `X_scVI` exists.
- After `CellANOVA`, check `denoised` plus `X_cellanova`.
- After `Concord`, check that a Concord object was returned and that `X_concord` exists.
- Before benchmarking, confirm that the batch and label columns both exist.
- After benchmarking, confirm that the results table is non-empty and includes the requested embedding keys.
- If only a representative smoke path was executed, state which backends and which benchmark subset were run empirically.

## Resource Map

- Use the source grounding notes for backend-specific outputs, return types, and prerequisite details.
- Use the backend selection notes when choosing among CPU, GPU, or control-pool branches.
- Use the notebook mapping notes to align notebook sections with the stage boundaries in this skill.
- Use the compatibility notes for dependency-sensitive and long-running backends.
