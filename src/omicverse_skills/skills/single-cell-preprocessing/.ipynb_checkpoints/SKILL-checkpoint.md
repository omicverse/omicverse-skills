---
name: omicverse-single-cell-preprocessing
description: Convert OmicVerse single-cell preprocessing and marker-discovery notebooks into a reusable, triggerable skill. Use when cleaning AnnData, choosing a preprocessing mode, building PCA and neighborhood graphs, clustering with Leiden, or extracting marker genes with OmicVerse.
---

# OmicVerse Single-Cell Preprocessing

## Goal

Turn raw or already-QC'd single-cell `AnnData` into a reusable execution spine: QC, preprocessing, PCA, graph construction, embeddings, Leiden clustering, and optional marker discovery. Keep marker discovery in the same skill because it only becomes useful after the same cluster-ready object exists.

## Quick Workflow

1. Inspect the input `AnnData` shape, layers, and whether a cluster key already exists.
2. If the data is raw, run `ov.pp.qc(...)` first; if it is already filtered, jump to preprocessing.
3. Choose one preprocessing branch with `ov.pp.preprocess(...)`, then keep the HVG output and build `scaled`, PCA, neighbors, and an embedding.
4. Run `ov.pp.leiden(...)` after the graph exists.
5. If you need marker genes, choose a `method` that matches the input representation, then call `ov.single.find_markers(...)` and `ov.single.get_markers(...)` or `ov.pl.markers_dotplot(...)`.
6. Validate the expected `obs`, `var`, `obsm`, `obsp`, and `uns` keys before writing output.

## Interface Summary

- `ov.pp.qc(adata, **kwargs)` handles QC and doublet filtering.
- `ov.pp.preprocess(adata, mode='shiftlog|pearson', target_sum=500000.0, n_HVGs=2000, organism='human', no_cc=False, batch_key=None, identify_robust=True)` normalizes and selects HVGs.
- `ov.pp.scale(adata, max_value=10, layers_add='scaled', to_sparse=True, **kwargs)` creates the `scaled` layer.
- `ov.pp.pca(adata, n_pcs=50, layer='scaled', inplace=True, **kwargs)` stores PCA in `scaled|original|X_pca`.
- `ov.pp.neighbors(adata, n_neighbors=15, n_pcs=None, use_rep=None, knn=True, random_state=0, n_jobs=None, method='umap', transformer=None, metric='euclidean', metric_kwds=mappingproxy({}), key_added=None, copy=False, **kwargs)` builds the graph.
- `ov.pp.umap(adata, **kwargs)`, `ov.pp.mde(adata, ...)`, `ov.pp.tsne(adata, ...)`, and `ov.pp.sude(adata, ...)` are alternative embedding branches.
- `ov.pp.leiden(adata, resolution=1.0, random_state=0, key_added='leiden', local_iterations=100, max_levels=10, device='cpu', symmetrize=None, **kwargs)` clusters the graph.
- `ov.single.find_markers(adata, groupby, method='cosg', n_genes=50, key_added=None, use_raw=None, layer=None, groups='all', reference='rest', corr_method='benjamini-hochberg', rankby_abs=False, tie_correct=False, pts=True, **kwargs)` finds markers.
- `ov.single.get_markers(adata, n_genes=10, key='rank_genes_groups', groups=None, return_type='dataframe', min_logfoldchange=None, min_score=None, min_pval_adj=None)` extracts marker tables.
- `ov.pl.markers_dotplot(...)` renders the marker summary after `find_markers`.

Read `references/source-grounding.md` before adding more interface-specific details.

## Stage Selection

- Use this branch selection to avoid replaying the whole notebook when the object is already partially prepared.
- Use the full pipeline when the input is raw counts and no graph exists yet.
- Skip directly to marker discovery when the object is already normalized, clustered, and only marker extraction or dotplotting is missing.
- Use `mode='shiftlog|pearson'` when you want the notebook's default preprocessing path.
- Use `method='wilcoxon'` for statistically tested markers on log-normalized data.
- Use `method='cosg'` when you want fast marker discovery from raw counts.
- Use `method='t-test'`, `method='t-test_overestim_var'`, or `method='logreg'` only when their assumptions match the data.

## Input Contract

- Start from `AnnData`.
- Preserve raw counts if you need a marker branch that expects them.
- Ensure the neighbor graph exists before `ov.pp.leiden(...)`.
- Ensure a cluster key such as `leiden` exists before `ov.single.find_markers(...)` and `ov.pl.markers_dotplot(...)`.
- Treat `ov.settings.cpu_gpu_mixed_init()` as a runtime toggle, not as a separate skill boundary.

## Minimal Execution Patterns

```python
import omicverse as ov

ov.settings.cpu_init()

adata = ov.pp.qc(adata, mode="mads", doublets=False, min_genes=0, min_cells=0)
adata = ov.pp.preprocess(
    adata,
    mode="shiftlog|pearson",
    n_HVGs=2000,
    target_sum=1e4,
    identify_robust=True,
)
adata.raw = adata
adata = adata[:, adata.var.highly_variable_features].copy()
ov.pp.scale(adata)
ov.pp.pca(adata, layer="scaled", n_pcs=50)
ov.pp.neighbors(adata, n_neighbors=15, n_pcs=50, use_rep="scaled|original|X_pca")
ov.pp.umap(adata)
ov.pp.leiden(adata, resolution=1.0)
ov.single.find_markers(adata, groupby="leiden", method="wilcoxon", key_added="rank_genes_groups")
markers = ov.single.get_markers(adata, key="rank_genes_groups", return_type="dict")
```

```python
ov.single.find_markers(adata, groupby="leiden", method="cosg", key_added="cosg_markers", pts=True)
ov.pl.markers_dotplot(adata, groupby="leiden", key="cosg_markers", n_genes=5, show=False)
```

## Constraints

- Do not hardcode notebook data paths or local environment names.
- Do not assume `scrublet` works on tiny fixtures; keep bounded smoke checks synthetic and lightweight.
- Keep `cosg` and the statistical methods separate: `cosg` expects raw counts, while `wilcoxon`, `t-test`, and `logreg` expect log-normalized data.
- Be explicit about `mode`, `method`, and `use_rep`; do not leave branch selection implicit.
- Keep the notebook's visualization polish cells optional unless the user explicitly wants the same figures.

## Validation

- Check `adata.layers["counts"]` after preprocessing.
- Check `adata.var["highly_variable_features"]` before PCA.
- Check `adata.obsm["scaled|original|X_pca"]` after PCA.
- Check `adata.obsm["X_umap"]` or another embedding key after the embedding step.
- Check `adata.obs["leiden"]` before marker discovery.
- Check `adata.uns["rank_genes_groups"]` or the selected marker key before `get_markers(...)` or `markers_dotplot(...)`.
- If you only validated a smoke path, say so; do not claim full-scale reproduction.

## Resource Map

- Read `references/branch-selection.md` when choosing QC, preprocessing, graph, embedding, clustering, or marker branches.
- Read `references/source-notebook-map.md` to trace notebook cells into this reusable skill.
- Read `references/source-grounding.md` for inspected signatures, source-level branch behavior, and compatibility notes.
- Read `references/compatibility.md` when notebook prose, parameter names, or current runtime behavior diverge.
