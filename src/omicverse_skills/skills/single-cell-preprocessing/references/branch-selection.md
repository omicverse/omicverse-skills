# Branch Selection

This notebook is one reusable skill because the later stages depend on the same `AnnData` state:

- QC produces a filtered input.
- preprocessing produces HVGs and normalized counts.
- PCA, neighbors, and embeddings produce the graph basis.
- Leiden adds cluster labels.
- marker discovery consumes the clustered object.

Do not split those stages into separate skills unless the user explicitly wants only marker extraction or only preprocessing.

## QC Branches

- Use `mode='seurat'` for threshold-based filtering.
- Use `mode='mads'` for robust outlier filtering.
- Use `doublets_method='scrublet'` when you want the notebook's original path and the fixture is large enough for scrublet to survive filtering.
- Use `doublets_method='sccomposite'` when scrublet is too brittle for the available data.
- Keep `filter_doublets=False` when you only need a QC callability smoke check.

## Preprocessing Branches

- Use `mode='shiftlog|pearson'` for the notebook's main path.
- Use `mode='shiftlog|seurat'` when you want Seurat-style HVG selection on log-normalized data.
- Use `mode='pearson|pearson'` or `mode='pearson|seurat'` when you want the Pearson-residual workflow.
- Keep `target_sum` and `n_HVGs` explicit when you need reproducible downstream `pca` and `neighbors` results.

## Graph And Embedding Branches

- Use `ov.pp.neighbors(..., method='umap')` for the default connectivity graph.
- Use `method='gauss'` only when you intentionally want Gaussian-kernel connectivities.
- Use `method='rapids'` only for the GPU backend path.
- Use `ov.pp.umap(...)` for the notebook's default visualization branch.
- Use `ov.pp.mde(...)`, `ov.pp.tsne(...)`, or `ov.pp.sude(...)` only when you need those alternative embeddings.
- Treat `ov.pl.mde(...)` as a plotting helper around `X_mde`, not a separate analysis stage.

## Clustering Branches

- Use `ov.pp.leiden(...)` after neighbors exist.
- The runtime backend is controlled by `ov.settings.mode`:
  - `cpu` uses the Scanpy branch.
  - `cpu-gpu-mixed` uses the mixed omicverse branch.
  - GPU mode routes to RAPIDS where available.
- Leave `resolution` explicit when you need comparable cluster granularity across runs.

## Marker Branches

- Use `method='wilcoxon'` when your data is log-normalized and you want a statistical test.
- Use `method='cosg'` when your data are raw counts and you want the notebook's faster marker branch.
- Use `method='t-test'`, `method='t-test_overestim_var'`, or `method='logreg'` only when their assumptions fit the question.
- Use `return_type='dict'` in `get_markers(...)` when you want a compact mapping for annotation or report generation.

## Why Keep One Skill

The notebook mixes multiple visibly independent plots, but they are all downstream of the same preprocessing spine. Splitting the plotting polish or marker export into a second skill would create a thin wrapper with the same input contract and would force callers to switch skills before the object is actually ready.
