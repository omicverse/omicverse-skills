---
name: omicverse-single-cell-clustering-backends
description: Run and compare OmicVerse single-cell clustering backends as a reusable, triggerable skill. Use when choosing between Leiden, Louvain, scICE, or GMM Gaussian mixture clustering on a prepared AnnData embedding, or when adapting a related OmicVerse clustering notebook into a repeatable workflow.
---

# OmicVerse Single-Cell Clustering Backends

## Goal

Turn the notebook's clustering comparison into one reusable job: choose and run a clustering backend on a preprocessed single-cell `AnnData`. Keep this skill focused on backend selection and output contracts after PCA and, when needed, graph construction.

## Quick Workflow

1. Inspect whether the input already has PCA or another embedding, and whether a neighbor graph already exists.
2. Choose a backend explicitly: `leiden`, `louvain`, `scICE`, or `GMM` for the notebook-covered paths.
3. Build a neighbor graph first for `leiden` and `louvain`.
4. Pass `use_rep` explicitly for embedding-based methods such as `GMM` and `scICE`.
5. Validate which `obs` columns were added, and treat ARI or figure rendering as optional downstream reporting rather than the core contract.

## Interface Summary

- `ov.utils.cluster(adata, method='leiden', use_rep='X_pca', random_state=1024, n_components=None, **kwargs)` dispatches to a clustering backend.
- Live source exposes more `method` values than the notebook covers: `leiden`, `louvain`, `kmeans`, `GMM`, `mclust`, `mclust_R`, `schist`, and `scICE`.
- `leiden` and `louvain` forward `**kwargs` to scanpy clustering on the current graph.
- `GMM` and `mclust` use the embedding named by `use_rep`, require `n_components`, and write both `mclust` and `gmm_cluster`.
- `scICE` returns a fitted model object and adds one or more `scICE_k*` columns to `adata.obs` through `add_to_adata(...)`.
- `scICE` runtime behavior depends on `resolution_range`, `n_steps`, `n_trials`, and `n_boot`.

## Stage Selection

- Use `leiden` for the default graph-based clustering path.
- Use `louvain` only when Louvain compatibility is specifically requested.
- Use `GMM` when the user wants soft-assignment-style Gaussian mixture clustering on an embedding instead of graph partitioning.
- Use `scICE` when the user wants consistency-guided cluster-count selection and is willing to pay the extra bootstrap cost.
- Keep the notebook's UMAP and ARI cells optional; they are reporting on top of clustering outputs, not the reusable skill boundary.

## Input Contract

- Start from an `AnnData` object that already has PCA or another usable embedding for embedding-based methods.
- Ensure a neighbor graph exists before `leiden` or `louvain`.
- Ensure the embedding named by `use_rep` exists before `GMM` or `scICE`.
- Pass `n_components` when using `GMM` or `mclust`.

## Minimal Execution Patterns

```python
import omicverse as ov

ov.pp.neighbors(adata, n_neighbors=15, n_pcs=50, use_rep="scaled|original|X_pca")
ov.utils.cluster(adata, method="leiden", resolution=1.0)
ov.utils.cluster(adata, method="louvain", resolution=1.0)
```

```python
ov.utils.cluster(
    adata,
    method="GMM",
    use_rep="scaled|original|X_pca",
    n_components=21,
)
```

```python
scice_model = ov.utils.cluster(
    adata,
    method="scICE",
    use_rep="scaled|original|X_pca",
    resolution_range=(4, 20),
    n_boot=50,
    n_steps=11,
)
```

## Constraints

- Do not use this skill as a preprocessing skill; it assumes the object is already prepared enough to cluster.
- Do not leave `method` implicit.
- Do not assume the notebook's cluster counts are stable defaults for other datasets.
- Treat `scICE_k*` outputs as data-dependent candidate solutions, not a fixed schema with one guaranteed suffix.
- Keep smoke and acceptance commands shell-agnostic.

## Validation

- Check that the intended `obs` column was added for the chosen backend.
- For `GMM`, check both `mclust` and `gmm_cluster`.
- For `scICE`, check that a model object was returned and that at least one `scICE_k*` column was added.
- If comparing backends, compute metrics such as ARI only after all requested label columns exist.
- If only a bounded smoke path was run, say which backends were executed and which were source-grounded only.

## Resource Map

- Use the branch selection notes when choosing a backend.
- Use the source grounding notes for current method branches and output columns.
- Use the notebook mapping notes to trace the notebook's sections into this reusable skill.
- Use the compatibility notes for dependency-sensitive or backend-sensitive paths.
