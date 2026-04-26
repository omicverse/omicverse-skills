---
name: omicverse-single-cell-monocle2-trajectory
description: Monocle2-style single-cell trajectory analysis on AnnData via the `ov.single.Monocle` class - DDRTree pseudotime + branch detection + per-gene differential test + BEAM branch-dependent gene discovery, plus the canonical pseudotime visualisations (`branch_streamplot`, `dynamic_heatmap`, `dynamic_trends`). Use when fitting a Monocle2 trajectory on an annotated AnnData, when deriving branch-aware gene trends with `dynamic_features`, or when reproducing `t_traj_monocle2_olsson`.
---

# OmicVerse Single-Cell — Monocle2 Trajectory

## Goal

Take a preprocessed annotated single-cell `AnnData` (typically myeloid / hematopoietic lineage) and fit a Monocle2-style trajectory: DDRTree-based ordering, pseudotime computation, branch-point detection, per-gene differential testing along pseudotime, and BEAM (Branch Expression Analysis Modeling) for branch-dependent genes. Visualise the result with the OmicVerse plotting stack (`branch_streamplot`, `dynamic_heatmap`, `dynamic_trends`) on top of the unified `dynamic_features` GAM-fitting backend.

`Monocle` here is a **wrapper** around the original Monocle2 logic re-implemented in Python; it operates on `AnnData` rather than `cell_data_set` and integrates with the rest of OmicVerse plotting / trajectory tooling. `dynamic_features` is the shared GAM backend used by every trajectory skill (Palantir, Slingshot, Monocle, etc.) so genes-along-pseudotime visualisations are interoperable.

## Quick Workflow

1. Load a preprocessed AnnData with cell-type labels (`obs[<celltype>]`) and counts in `.X`. The Olsson 2016 hematopoietic dataset is the canonical demo.
2. Instantiate: `mono = Monocle(adata)`. Stores a copy internally; `mono.adata` is the working object.
3. **Preprocess + ordering genes**: `mono.preprocess()` (variance-stabilises + filters) → `mono.select_ordering_genes(max_genes=1000)` (high-variance + DEG-driven gene selection); inspect with `mono.plot_ordering_genes(figsize=(5, 4))`. Result: `mono.adata.var['use_for_ordering']` boolean mask.
4. **Fit trajectory**: by default `mono.fit_trajectory()` runs DDRTree on the ordering-gene subspace and writes `obs['Pseudotime']`, `obs['State']` (DDRTree-derived branch state). Call directly OR run via `mono.plot_trajectory(...)` which auto-fits on first call.
5. **Plot**: `mono.plot_trajectory(color_by='State', cell_size=2.4, cell_link_size=0.5, show_branch_points=False, figsize=(6, 4))` for the canonical Monocle2 tree-on-DDRTree-2D figure.
6. **Branch streamplot** (alternative pseudotime view): `ov.pl.branch_streamplot(mono.adata, group_key=<celltype>, pseudotime_key='Pseudotime')` — colour-coded streamlines along pseudotime.
7. **Differential gene test along pseudotime**: re-instantiate on the ordering-gene subset for speed: `mono_ord = Monocle(mono.adata[:, ordering_genes].copy())` then `de = mono_ord.differential_gene_test(cores=-1)`. Filter `de[(de['qval']<0.01) & (de['status']=='OK')]` for significant trends.
8. **Heatmap of top trending genes**: `ov.pl.dynamic_heatmap(mono.adata, pseudotime='Pseudotime', var_names=top40, cell_annotation='State', use_fitted=True, cell_bins=200, ...)`.
9. **Per-gene trends**: `res = ov.single.dynamic_features(mono.adata, genes=marker_genes, pseudotime='Pseudotime', store_raw=True, raw_obs_keys=['State'])` → `ov.pl.dynamic_trends(res, genes=marker_genes, add_point=True, point_color_by='State', ...)`.
10. **Branch-aware trends** (compare two trajectories that share an early trunk): pass `groupby=<celltype>, groups=['Gmp', 'LK']` to `dynamic_features`, then plot with `compare_groups=True, split_time=<trunk-end-time>, shared_trunk=True`.
11. **BEAM branch-dependent genes**: `mono_ord.BEAM(branch_point=1, cores=-1) → DataFrame` with `pval`, `qval` per gene; combine with `dynamic_features(groupby=..., groups=branch_subtypes)` + `dynamic_trends(compare_groups=True, ...)` to visualise the branching pattern.

## Interface Summary

```python
ov.single.Monocle(adata: AnnData)
# Stores adata in `self.adata` (a copy). Subsequent methods mutate self.adata.
```

Class methods (chained or independent):
- `mono.preprocess()` — variance-stabilising normalisation; canonical Monocle2 prep.
- `mono.select_ordering_genes(max_genes=1000)` — populates `var['use_for_ordering']`. The DEG-driven side combines with high-variance genes.
- `mono.plot_ordering_genes(figsize=...)` — diagnostic plot of the ordering-gene selection.
- `mono.fit_trajectory(...)` — DDRTree backbone fit; writes `obs['Pseudotime']`, `obs['State']`. Triggered automatically by `plot_trajectory` on first call.
- `mono.plot_trajectory(color_by='State', cell_size=..., cell_link_size=..., show_branch_points=..., figsize=..., palette=..., ...)` — Monocle2's signature DDRTree-2D scatter.
- `mono.differential_gene_test(cores=-1, ...) → pd.DataFrame` — per-gene VGAM-style tests for variation along pseudotime; `cores=-1` uses all cores. Output columns: `pval`, `qval`, `status` (`'OK'` / `'FAIL'`), `gene_id`.
- `mono.BEAM(branch_point=1, cores=-1, ...) → pd.DataFrame` — Branch Expression Analysis Modeling: per-gene tests for branch-dependent expression at the named branch point. Same output schema as `differential_gene_test`.

Module-level (not method on `Monocle`):
- `ov.single.dynamic_features(adata|dict[name→adata], genes, pseudotime='pseudotime', *, groupby=None, groups=None, layer=None, use_raw=False, subsets=None, weights=None, distribution='normal', link='identity', n_splines=8, spline_order=3, grid_size=200, confidence_level=0.95, min_cells=20, min_variance=1e-08, store_raw=False, raw_obs_keys=None, key_added='dynamic_features', verbose=True) → DynamicFeaturesResult`. Generic GAM backend usable across any pseudotime field.
- `ov.pl.branch_streamplot(adata, *, group_key, pseudotime_key, show=True, ...)`.
- `ov.pl.dynamic_heatmap(adata, *, pseudotime, var_names, cell_annotation=None, use_cell_columns=False, use_fitted=True, cell_bins=200, smooth_window=..., fitted_window=..., figsize=..., show_row_names=True, standard_scale='var', cmap='RdBu_r', order_by='peak', show=True, ...)`.
- `ov.pl.dynamic_trends(res, *, genes, add_point=True, point_color_by=None, line_style_by=None, compare_features=False, compare_groups=False, split_time=None, shared_trunk=True, figsize=..., legend_loc=..., ncols=..., title=..., ...)`.

## Boundary

**Inside scope:**
- Full Monocle2 workflow on a preprocessed annotated AnnData.
- DEG along pseudotime (`differential_gene_test`).
- BEAM branch-dependent gene discovery.
- Pseudotime-along-genes visualisation (`dynamic_features` + `dynamic_trends` + `dynamic_heatmap`).
- Branch streamplot (`ov.pl.branch_streamplot`).

**Outside scope — separate skill:**
- Diffusion / Slingshot / Palantir trajectory inference — see `omicverse-single-cell-trajectory-inference` (`TrajInfer` wrapper).
- VIA trajectory inference (with or without RNA velocity) — see `omicverse-single-cell-via-trajectory`.
- scTour pseudotime — see `omicverse-single-cell-sctour-trajectory` (separate existing skill).
- CellRank-style fate maps from RNA velocity — see `omicverse-single-cell-cellrank-fate`.
- Trajectory inference on multi-modal (e.g. ATAC) data — out of scope; use `dynamic_features(layer='atac')` only after fitting a separate trajectory.
- Preprocessing / clustering / annotation — separate existing skills.

## Branch Selection

**`max_genes` for `select_ordering_genes`**
- 1000 (tutorial default) — good balance of biological coverage and DDRTree compute time on cohorts of 1k–10k cells.
- 500 — for small cohorts (<2k cells); too many ordering genes cause DDRTree to over-branch.
- 2000 — for large cohorts (>20k cells) with rich biology; check that `differential_gene_test` runtime stays acceptable.

**`differential_gene_test` cores**
- `cores=-1` uses all cores. Per-gene independent tests parallelise well; on 1000 genes × 8 cores expect ~30 s.
- Limit to physical cores (not hyperthreaded) for best throughput.

**`BEAM` branch_point**
- `branch_point=1` — first branch in the DDRTree (tutorial default for the Olsson hematopoietic dataset). Inspect `mono.adata.uns` or the `plot_trajectory` figure to identify the right branch number.
- For multi-branch trajectories, run BEAM per branch point separately and concatenate; never run on `branch_point=0` (the trunk) — that's not a branch.

**`dynamic_features` `distribution` / `link`**
- `'normal' / 'identity'` (default) — the right choice for log-normalised expression (post `sc.pp.log1p`). Tutorial default.
- `'gamma' / 'log'` — for raw counts directly; rarely needed since most pipelines log-normalise first.
- `'binomial' / 'logit'` — only for binary on/off detection traces, not log-CPM expression.

**`dynamic_features` `n_splines` / `spline_order`**
- 8 splines, order 3 (cubic) — tutorial default; smooth without over-fitting.
- Increase to 12–16 splines for traces with multiple inflection points (e.g. transient peaks).
- Decrease to 4–6 for cohorts with sparse pseudotime coverage.

**`dynamic_trends` modes**
- Default: one panel per gene, with cells coloured by `point_color_by`.
- `compare_features=True` — overlay multiple gene curves on the same axis; `line_style_by='features'` differentiates them. Useful for showing ordered marker progression.
- `compare_groups=True` — branch-aware: requires `groupby='subtype'` (or similar) at `dynamic_features` time and `groups=[branch_a, branch_b]`. Pass `split_time=<trunk-to-branch transition>` and `shared_trunk=True` to draw a single trunk that splits.

**`dynamic_heatmap` `use_fitted=True` vs raw smoothing**
- `use_fitted=True` — uses GAM-fitted curves (smoothest; cleanest for publication).
- `use_fitted=False` — uses raw smoothed bins; better for spotting outlier cells.
- `cell_bins=200` is default; reduce to 50–100 for small cohorts (<1k cells).

**Subsetting before `differential_gene_test`**
- Re-instantiate `Monocle(mono.adata[:, ordering_genes].copy())` before the test — the test runs over all genes by default, and on a 20k-gene matrix with 8 cores it's ~minutes. Subsetting to ~1000 ordering genes drops it to <1 minute.

## Input Contract

- `AnnData` with raw counts (or normalised log-counts) in `.X`. `obs` should carry the cell-type column you'll pass to `branch_streamplot` (`group_key`) and trend plots (`point_color_by`).
- For Olsson-style data: `obs` includes a `subtype` column (`'Gmp'`, `'LK'`, `'Cmp'`, etc.); the WT subset is the canonical demo (`adata_wt = adata[adata.obs['genotype'] == 'WT']`).
- `var` has gene symbols (or IDs); `differential_gene_test` indexes by `var_names`.
- `obs['Pseudotime']` and `obs['State']` are populated by `fit_trajectory`. Don't pass these in; they'll be overwritten.

## Minimal Execution Patterns

```python
import anndata as ad
import omicverse as ov
from omicverse.single import Monocle
import matplotlib.pyplot as plt

ov.plot_set(font_path='Arial')

# 1) Load preprocessed annotated AnnData (Olsson WT subset)
adata_wt = ad.read_h5ad('olsson_wt.h5ad')

# 2) Fit Monocle2
mono = Monocle(adata_wt)
mono.preprocess()
mono.select_ordering_genes(max_genes=1000)
mono.plot_ordering_genes(figsize=(5, 4)); plt.show()

# 3) Plot trajectory (auto-fits on first call)
mono.plot_trajectory(
    color_by='State', cell_size=2.4, cell_link_size=0.5,
    show_branch_points=False, figsize=(6, 4),
); plt.show()

# 4) Branch streamplot view
fig, ax = ov.pl.branch_streamplot(
    mono.adata, group_key='subtype', pseudotime_key='Pseudotime', show=False,
); plt.show()

# 5) DEG along pseudotime — re-instantiate on ordering subset for speed
ordering = mono.adata.var_names[mono.adata.var['use_for_ordering']].tolist()
mono_ord = Monocle(mono.adata[:, ordering].copy())
de = mono_ord.differential_gene_test(cores=-1)
sig = de[(de['qval'] < 0.01) & (de['status'] == 'OK')]
top40 = sig.sort_values('pval').head(40).index.tolist()

# 6) Heatmap of top trending genes (GAM-fitted)
g = ov.pl.dynamic_heatmap(
    mono.adata,
    pseudotime='Pseudotime',
    var_names=top40,
    cell_annotation='State',
    use_fitted=True,
    cell_bins=200,
    figsize=(7, 7),
    show_row_names=True,
    standard_scale='var', cmap='RdBu_r', order_by='peak',
    show=False,
); plt.show()

# 7) Per-gene trends with cell points
markers = [g for g in ['Gfi1', 'Irf8', 'Elane', 'Prtn3', 'Mpo', 'Car2']
           if g in mono.adata.var_names]
res = ov.single.dynamic_features(
    mono.adata, genes=markers, pseudotime='Pseudotime',
    store_raw=True, raw_obs_keys=['State'],
)
ov.pl.dynamic_trends(
    res, genes=markers,
    add_point=True, point_color_by='State',
    figsize=(4, 3.5), legend_loc='right margin', legend_fontsize=8,
); plt.show()
```

```python
# Branch-aware trends — compare two branches that share an early trunk
import numpy as np

branch_subtypes = ['Gmp', 'LK']
res_branch = ov.single.dynamic_features(
    mono.adata,
    genes=['Gfi1', 'Irf8', 'Elane', 'Car2'],
    pseudotime='Pseudotime',
    groupby='subtype',
    groups=branch_subtypes,
    store_raw=True,
)
trunk_mask = mono.adata.obs['subtype'].astype(str).isin(['Cmp'])
split_time = float(np.nanmedian(mono.adata.obs.loc[trunk_mask, 'Pseudotime']))
ov.pl.dynamic_trends(
    res_branch, genes=['Gfi1', 'Irf8', 'Elane', 'Car2'],
    compare_groups=True,
    split_time=split_time, shared_trunk=True,
    add_point=True, point_color_by='group',
    figsize=(4.6, 3), ncols=2,
)

# BEAM at branch point 1
beam = mono_ord.BEAM(branch_point=1, cores=-1)
sig_beam = beam[beam['qval'] < 0.01].sort_values('qval')
print(f'BEAM hits: {len(sig_beam)} / {len(beam)}')
top_branch = sig_beam.head(4).index.tolist()
```

## Validation

- After `select_ordering_genes`: `mono.adata.var['use_for_ordering'].sum()` should be close to `max_genes` (slight under is OK if some genes failed variance filter; large under means filter was too aggressive).
- After `fit_trajectory` (or `plot_trajectory`): `mono.adata.obs['Pseudotime']` and `mono.adata.obs['State']` are populated, no NaNs.
- `differential_gene_test`: `(de['status'] == 'OK').sum() / len(de)` should be > 0.9. Many `'FAIL'` rows indicate genes with insufficient signal after filtering — usually fine, just lowers the universe of testable genes.
- `(de['qval'] < 0.01).sum()` typically yields 50–500 hits on Olsson-scale cohorts; far fewer indicates the trajectory is poorly defined.
- BEAM hits should be a *subset* of `differential_gene_test` hits (with the additional branch-dependence constraint). If BEAM returns *more* hits than DEG, something is mis-specified (likely the wrong `branch_point`).
- `dynamic_features` always returns a `DynamicFeaturesResult` even when some genes failed; the `.failed` attribute lists failures with reasons (insufficient cells, low variance, GAM fit issues).
- `dynamic_heatmap` over `top40` genes: rows should show clear pseudotime gradients (the whole point of selecting on q<0.01); if the heatmap looks random, the trajectory or pseudotime is degenerate.

## Resource Map

- See [`reference.md`](reference.md) for compact copy-paste snippets.
- See [`references/source-grounding.md`](references/source-grounding.md) for verified `Monocle` / `dynamic_features` / `branch_streamplot` / `dynamic_heatmap` / `dynamic_trends` signatures.
- For other trajectory backends (diffusion / slingshot / palantir), see `omicverse-single-cell-trajectory-inference`; for VIA, see `omicverse-single-cell-via-trajectory`; for scTour, see `omicverse-single-cell-sctour-trajectory`.
- For RNA velocity-based fate maps, see `omicverse-single-cell-cellrank-fate`.

## Examples
- "Fit Monocle2 on the Olsson WT subset, plot the DDRTree trajectory by State, and run differential_gene_test on the ordering-gene subset."
- "Run BEAM at branch point 1 and visualise the top-4 branch-dependent genes with `compare_groups=True, shared_trunk=True`."
- "Build a 40-gene dynamic heatmap of top pseudotime-significant genes ordered by peak, GAM-fitted."
- "Show per-gene trends for Gfi1 / Irf8 / Elane along Pseudotime with cell points coloured by State."

## References
- Tutorial notebook: [`t_traj_monocle2_olsson.ipynb`](https://omicverse.readthedocs.io/en/latest/Tutorials-single/t_traj_monocle2_olsson/) — Olsson 2016 hematopoietic walkthrough.
- Live API verified — see [`references/source-grounding.md`](references/source-grounding.md).
