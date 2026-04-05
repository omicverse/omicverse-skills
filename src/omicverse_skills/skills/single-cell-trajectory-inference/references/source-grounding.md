# Source Grounding

This skill is grounded in the current local OmicVerse source tree and live interface inspection against the importable checkout.

## Primary Source Files

- `omicverse/single/_traj.py`
- `omicverse/single/_pyslingshot.py`
- `omicverse/utils/_paga.py`
- `omicverse/external/palantir/core.py`
- `omicverse/external/palantir/presults.py`

## Inspected Interfaces

### `ov.single.TrajInfer`

Signature:

```text
(adata: anndata.AnnData, basis: str = 'X_umap', use_rep: str = 'X_pca', n_comps: int = 50, n_neighbors: int = 15, groupby: str = 'clusters')
```

Source behavior:

- `basis` is the embedding key used for visualization-facing methods
- `use_rep` is the representation key consumed by trajectory computation
- `groupby` is the cluster key used to resolve origin and terminal labels

### `TrajInfer.inference`

Signature:

```text
(self, method: str = 'palantir', **kwargs)
```

Observed source branches in `omicverse/single/_traj.py`:

- `method='palantir'`
- `method='diffusion_map'`
- `method='slingshot'`
- `method='sctour'`

Branch behavior:

- `diffusion_map` writes `dpt_pseudotime`
- `slingshot` writes `slingshot_pseudotime`
- `palantir` writes `palantir_pseudotime` and returns the Palantir result object
- `sctour` writes `sctour_pseudotime`, `X_TNODE`, and `X_VF`

### `Slingshot.fit`

Source signature from `omicverse/single/_pyslingshot.py`:

```text
fit(self, num_epochs=10, debug_axes=None)
```

Observed behavior:

- `num_epochs` controls repeated curve-fitting epochs
- `debug_axes` only drives the optional visualization path
- the implementation imports `pcurvepy2.PrincipalCurve` during curve construction

### `ov.utils.cal_paga`

Signature:

```text
(adata, groups=None, vkey='velocity', use_time_prior=True, root_key=None, end_key=None, threshold_root_end_prior=None, minimum_spanning_tree=True, copy=False)
```

Observed source behavior:

- requires `adata.uns['neighbors']`
- auto-detects `groups` as `clusters` or `louvain` only when left unset
- writes `adata.uns['paga']`
- uses `use_time_prior` to bias transitions with a pseudotime key when supplied

### `ov.utils.plot_paga`

Signature:

```text
(adata, basis=None, vkey='velocity', color=None, layer=None, title=None, threshold=None, layout=None, layout_kwds=None, init_pos=None, root=0, labels=None, single_component=False, dashed_edges='connectivities', solid_edges='transitions_confidence', transitions='transitions_confidence', node_size_scale=1, node_size_power=0.5, edge_width_scale=0.4, min_edge_width=None, max_edge_width=2, arrowsize=15, random_state=0, pos=None, node_colors=None, normalize_to_color=False, cmap=None, cax=None, cb_kwds=None, add_pos=True, export_to_gexf=False, plot=True, use_raw=None, size=None, groups=None, components=None, figsize=None, dpi=None, show=None, save=None, ax=None, ncols=None, scatter_flag=None, **kwargs)
```

Notebook-relevant parameters:

- `basis`
- `color`
- `title`
- `min_edge_width`
- `node_size_scale`
- `show`

### `run_palantir`

Signature:

```text
(data, early_cell, terminal_states=None, knn=30, num_waypoints=1200, n_jobs=-1, scale_components=True, use_early_cell_as_start=False, max_iterations=25, eigvec_key='DM_EigenVectors_multiscaled', pseudo_time_key='palantir_pseudotime', entropy_key='palantir_entropy', fate_prob_key='palantir_fate_probabilities', save_as_df=None, waypoints_key='palantir_waypoints', seed=20)
```

Notebook-relevant parameters:

- `terminal_states`
- `knn`
- `num_waypoints`
- `n_jobs`

### `select_branch_cells`

Signature:

```text
(ad, pseudo_time_key='palantir_pseudotime', fate_prob_key='palantir_fate_probabilities', q=0.01, eps=0.01, masks_key='branch_masks', save_as_df=None)
```

Notebook-relevant parameters:

- `eps`
- `masks_key`

### `compute_gene_trends`

Signature:

```text
(ad, lineages=None, masks_key='branch_masks', expression_key=None, pseudo_time_key='palantir_pseudotime', gene_trend_key='gene_trends', save_as_df=None, **kwargs)
```

Notebook-relevant parameters:

- `lineages`
- `masks_key`
- `expression_key`
- `gene_trend_key`

## Runtime Findings That Affect This Skill

- The Slingshot branch imports `pcurvepy2` and fails without it.
- The Palantir gene-trend branch imports `mellon` and fails without it.
- The Palantir wrapper currently calls `run_magic_imputation(self.adata)` without surfacing the helper's `n_jobs` argument, so constrained local smoke harnesses may need a single-process patch to validate the branch safely.

## Empirical Smoke Coverage Used For This Skill

- synthetic `diffusion_map` run to confirm `dpt_pseudotime` plus PAGA output
- synthetic Palantir run with bounded `knn`, `num_waypoints`, and a single-process MAGIC patch to confirm `palantir_pseudotime`, `palantir_entropy`, `palantir_fate_probabilities`, `branch_masks`, and PAGA output

The Palantir gene-trend branch and Slingshot branch were source-grounded but not fully reproduced in the smoke harness because the current environment lacks the needed optional dependencies.
