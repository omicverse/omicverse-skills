---
name: omicverse-single-cell-cellrank-fate
description: CellRank fate maps from RNA velocity. Combine VelocityKernel + ConnectivityKernel into a transition matrix, fit a GPCCA estimator, predict terminal states, and produce per-cell fate probabilities. Visualise with `ov.pl.branch_streamplot` and feed branch-resolved gene-trends into `ov.single.dynamic_features` / `ov.pl.dynamic_trends` / `ov.pl.dynamic_heatmap`. Use after RNA velocity is computed (scvelo / dynamo / latentvelo / graphvelo) and before reporting fate probabilities or marker dynamics.
---

# OmicVerse Single-Cell — CellRank Fate Maps from RNA Velocity

## Goal

Take a single-cell `AnnData` that already has RNA velocity computed (`adata.layers['velocity']`, `adata.uns['velocity_graph']`, etc.) and produce a fate map: terminal-state assignments, per-cell fate probabilities, and branch-aware marker dynamics. CellRank does the velocity-driven random-walk side; OmicVerse plotting + `dynamic_features` does the visualisation and gene-trend side.

This skill is the **post-velocity** skill: assumes velocity is already computed (use `single-cell-rna-velocity` for that). It pairs CellRank's VelocityKernel + ConnectivityKernel mix with the same `dynamic_features` GAM backend used by Monocle / VIA / etc., so plot helpers (`branch_streamplot`, `dynamic_trends`, `dynamic_heatmap`) are interoperable across all trajectory skills.

## Quick Workflow

1. Confirm velocity is already in `adata.layers['velocity']` and `adata.uns['velocity_graph']` is populated. The pancreas endocrinogenesis dataset (`scv.datasets.pancreas`) is the canonical demo.
2. Build the **transition matrix** by mixing kernels:
   ```python
   vk = cr.kernels.VelocityKernel(adata).compute_transition_matrix()
   ck = cr.kernels.ConnectivityKernel(adata).compute_transition_matrix()
   kernel = 0.8 * vk + 0.2 * ck
   ```
   Velocity gives directionality; connectivity smooths over noise. The 80/20 mix is canonical.
3. Fit a **GPCCA estimator**: `g = cr.estimators.GPCCA(kernel)`. `g.compute_macrostates(n_states=6, cluster_key='clusters')` clusters the spectral decomposition; the resulting macrostates are candidates for terminal states.
4. **Predict terminal states**: `g.predict_terminal_states(n_states=4)`. Fewer terminals = more confidence per terminal; pick `n_states` from biological knowledge (4 for pancreas endocrine: Alpha / Beta / Delta / Epsilon, with the trunk macrostates excluded).
5. **Fate probabilities**: `g.fate_probabilities` is an `(n_cells, n_terminals)` DataFrame; pull a per-terminal probability and write to `adata.obs` for plotting:
   `adata.obs['beta_fate'] = np.asarray(g.fate_probabilities['Beta']).ravel()`.
6. **Branch streamplot** of pseudotime by cluster: `ov.pl.branch_streamplot(adata, group_key='clusters', pseudotime_key='development_pseudotime', trunk_groups=[...], branch_center=0.62, ...)`.
7. **Embedding-coloured fate**: `ov.pl.embedding(adata, basis='X_umap', color=['clusters', 'beta_fate'], cmap='Reds', frameon='small')`.
8. **Per-fate gene trends**: filter cells to the fate of interest (`adata_beta = adata[adata.obs['beta_fate'] > 0.15]`), then run `ov.single.dynamic_features(adata_beta, genes=[...], pseudotime='development_pseudotime', layer='Ms', ...)` and plot with `ov.pl.dynamic_trends(compare_features=True, ...)`. The `Ms` layer (smoothed-spliced from scvelo) is the appropriate input for velocity-derived workflows.
9. **Branch-aware comparison** between two terminal-state fates (e.g. Alpha vs Beta in pancreas): pass `groupby='clusters', groups=branch_clusters` to `dynamic_features`, then plot with `compare_groups=True, split_time=<trunk-end>, shared_trunk=True`.
10. **Multi-module dynamic heatmap**: organise marker genes into named programs (`{'Endocrine progenitor': [...], 'Alpha fate': [...], ...}`) and pass to `ov.pl.dynamic_heatmap(var_names=<dict>, ...)` for a panel-style figure.

## Interface Summary

External (CellRank — `pip install cellrank`):
- `cr.kernels.VelocityKernel(adata).compute_transition_matrix()`.
- `cr.kernels.ConnectivityKernel(adata).compute_transition_matrix()`.
- Kernel arithmetic: `0.8 * vk + 0.2 * ck` returns a combined kernel (CellRank operator-overloads).
- `cr.estimators.GPCCA(kernel)`.
  - `.compute_macrostates(n_states=N, cluster_key='clusters')`.
  - `.predict_terminal_states(n_states=k)`.
  - `.fate_probabilities` (DataFrame; `(n_cells, n_terminals)`).
  - `.macrostates`, `.terminal_states` — categorical columns.

OmicVerse-side (this skill's contribution):
- `ov.pl.branch_streamplot(adata, *, group_key, pseudotime_key, trunk_groups=[...], branch_center=0.5, figsize=..., xlabel=..., show=True)` — stream + branch view of per-cluster pseudotime.
- `ov.pl.embedding(adata, *, basis='X_umap', color=[...], cmap='Reds', frameon='small')` — generic embedding plot; colour by fate probability.
- `ov.single.dynamic_features(adata, genes, pseudotime, *, layer='Ms', groupby=None, groups=None, n_splines=8, store_raw=True, raw_obs_keys=...)` — same API as the trajectory skills, but `layer='Ms'` is the key choice for velocity workflows.
- `ov.pl.dynamic_trends(res, *, genes, compare_features=False, compare_groups=False, split_time=None, shared_trunk=True, add_point=True, point_color_by=..., line_style_by=..., ...)`.
- `ov.pl.dynamic_heatmap(adata, *, var_names: dict|list, pseudotime, layer='Ms', cell_annotation, ..., use_fitted=True, ...)`.

## Boundary

**Inside scope:**
- Building VelocityKernel + ConnectivityKernel mixes; choosing weights.
- GPCCA macrostate computation + terminal-state prediction.
- Pulling fate probabilities into `adata.obs` and plotting.
- Branch-aware gene-trend plots filtered by fate probability or terminal cluster.
- Multi-module dynamic heatmap.

**Outside scope — separate skill:**
- Computing RNA velocity itself — see `single-cell-rna-velocity` (existing skill; covers scvelo / dynamo / latentvelo / graphvelo).
- Other estimators: CFLARE, Schur — only GPCCA is documented here (the canonical CellRank 2 estimator).
- VIA / Monocle / diffusion / slingshot / palantir / scTour — see the per-method skills.
- Building the embedding (UMAP / FA2 / PHATE) — preprocessing skill.

## Branch Selection

**Kernel mix (`0.8 * vk + 0.2 * ck`)**
- Default tutorial mix: 80 % velocity, 20 % connectivity. Velocity drives directionality; connectivity bridges sparse-graph regions.
- 100 % velocity (`vk` alone): pure velocity-driven; works on dense, well-estimated velocities (e.g. Smart-seq2). Often too noisy on shallow 10x data.
- 50/50 mix: when velocity confidence is moderate; gives connectivity more weight to smooth out noise.
- 100 % connectivity (`ck` alone): no velocity at all — equivalent to a graph-based pseudotime; use only as a sanity baseline.

**`n_states` for `compute_macrostates`**
- Run with `n_states=6–10` first; inspect `g.macrostates.cat.categories`. Macrostates that look like trunks (high in-flow, low out-flow) are *not* terminal candidates.
- Then `predict_terminal_states(n_states=k)` where `k` is your biological estimate (e.g. 4 endocrine fates for pancreas).

**`n_states` for `predict_terminal_states`**
- Use biological knowledge: count distinct terminal cell types you expect from the cohort.
- Pancreas endocrinogenesis: 4 (Alpha / Beta / Delta / Epsilon) typically.
- Hematopoietic: 4–6 depending on cohort scope.
- Wrong: choosing `n_states` based on macrostate count alone — macrostates include trunks.

**Fate probability cutoff for filtering**
- 0.15 (tutorial default) — keep cells with at least 15 % probability of the chosen fate.
- 0.05 — looser; includes more cells (better for noisy GAM fits).
- 0.5 — strict; only confident-fate cells (cleaner trends but fewer cells, may under-fit).

**`pseudotime_key`** for `branch_streamplot` / `dynamic_features`
- `'development_pseudotime'` — when scvelo wrote it (`scv.tl.velocity_pseudotime`).
- `'velocity_pseudotime'` — same data, different naming convention.
- `'dpt_pseudotime'` — diffusion pseudotime (alternative; produces a different ordering).
- For CellRank-driven analysis, prefer `velocity_pseudotime` so the timing axis is consistent with the kernel.

**`layer='Ms'`** is critical: the smoothed spliced layer from scvelo. Don't pass `'spliced'` (raw, noisy) or `'velocity'` (vector, not expression). Without `layer='Ms'`, `dynamic_features` falls back to `adata.X` which is usually the wrong scale for post-velocity plotting.

**Branch-aware vs single-fate trends**
- Single-fate (`adata_beta = adata[fate > threshold]`): clean trends per terminal; doesn't show divergence.
- Branch-aware (`groupby='clusters', groups=[branch_a, branch_b]`, `compare_groups=True, shared_trunk=True`): shows where two fates diverge; needs a sensible `split_time` (median pseudotime of the trunk clusters).

## Input Contract

- `AnnData` with RNA velocity already computed:
  - `adata.layers['velocity']` populated.
  - `adata.uns['velocity_graph']` (or equivalent) for VelocityKernel.
  - `adata.obs['development_pseudotime']` (or `velocity_pseudotime`) for branch_streamplot / dynamic_features.
  - `adata.layers['Ms']` (smoothed spliced) for `dynamic_features(layer='Ms', ...)`.
  - `adata.obs['clusters']` (or any cell-type column) populated.
- `pip install cellrank` (the omicverse install does not ship CellRank by default).
- Embedding (`X_umap` / `X_pca`) for the visualisation steps.

## Minimal Execution Patterns

```python
import omicverse as ov
import scanpy as sc
import scvelo as scv
import cellrank as cr
import numpy as np
import matplotlib.pyplot as plt

ov.plot_set(font_path='Arial')

# Assume velocity has already been computed (see single-cell-rna-velocity skill)
# adata.layers['velocity'], adata.uns['velocity_graph'],
# adata.obs['development_pseudotime'], adata.layers['Ms'], adata.obs['clusters']

# 1) Build the kernel mix
vk = cr.kernels.VelocityKernel(adata).compute_transition_matrix()
ck = cr.kernels.ConnectivityKernel(adata).compute_transition_matrix()
kernel = 0.8 * vk + 0.2 * ck

# 2) GPCCA macrostates -> terminal states
g = cr.estimators.GPCCA(kernel)
g.compute_macrostates(n_states=6, cluster_key='clusters')
print('macrostates:', list(g.macrostates.cat.categories))

g.predict_terminal_states(n_states=4)
print('terminals:', list(g.terminal_states.cat.categories))

# 3) Per-cell fate probabilities -> obs column
adata.obs['beta_fate'] = np.asarray(g.fate_probabilities['Beta']).ravel()
ov.pl.embedding(adata, basis='X_umap',
                color=['clusters', 'beta_fate'],
                cmap='Reds', frameon='small')

# 4) Branch streamplot
fig, ax = ov.pl.branch_streamplot(
    adata, group_key='clusters',
    pseudotime_key='development_pseudotime',
    trunk_groups=['Ductal', 'Ngn3 low EP', 'Ngn3 high EP', 'Pre-endocrine'],
    branch_center=0.62,
    figsize=(6, 4),
    xlabel='Development pseudotime',
    show=False,
); plt.show()
```

```python
# Per-fate gene trends (Beta lineage)
adata_beta = adata[adata.obs['beta_fate'] > 0.15].copy()

beta_dyn = ov.single.dynamic_features(
    adata_beta,
    genes=['Pdx1', 'Ins1', 'Ins2', 'Sox9'],
    pseudotime='development_pseudotime',
    layer='Ms',
    distribution='normal', link='identity',
    n_splines=8, store_raw=True,
    raw_obs_keys=['clusters'],
)
ov.pl.dynamic_trends(
    beta_dyn, genes=['Pdx1', 'Ins1', 'Ins2', 'Sox9'],
    compare_features=True,
    add_point=True, point_color_by='clusters',
    line_style_by='features',
    figsize=(6.2, 3.2), linewidth=2.2,
)
```

```python
# Branch-aware comparison (Alpha vs Beta on the shared endocrine trunk)
branch_clusters = ['Alpha', 'Beta']
split_mask = adata.obs['clusters'].astype(str).isin(['Ngn3 high EP', 'Pre-endocrine'])

cr_branch_dyn = ov.single.dynamic_features(
    adata,
    genes=['Gcg', 'Ins1', 'Ins2', 'Pdx1'],
    pseudotime='development_pseudotime',
    layer='Ms',
    groupby='clusters',
    groups=branch_clusters,
    n_splines=8, store_raw=True,
)
split_time = float(np.nanmedian(
    adata.obs.loc[split_mask, 'development_pseudotime']
))
ov.pl.dynamic_trends(
    cr_branch_dyn, genes=['Gcg', 'Ins2', 'Pdx1'],
    compare_groups=True,
    split_time=split_time, shared_trunk=True,
    add_point=True, point_color_by='group',
    figsize=(4.2, 3), linewidth=2.2, ncols=3,
    title='CellRank branch-aware marker trends',
)
```

```python
# Multi-module dynamic heatmap of marker programs
modules = {
    'Endocrine progenitor': ['Sox9', 'Neurog3', 'Fev'],
    'Alpha fate':           ['Gcg', 'Arx'],
    'Beta fate':            ['Pax4', 'Ins1', 'Ins2', 'Pdx1'],
    'Delta fate':           ['Sst', 'Hhex'],
}
modules = {k: [g for g in v if g in adata.var_names] for k, v in modules.items()}
modules = {k: v for k, v in modules.items() if v}

ov.pl.dynamic_heatmap(
    adata,
    var_names=modules,
    pseudotime='development_pseudotime',
    layer='Ms',
    cell_annotation='clusters',
    cell_bins=180, smooth_window=17, fitted_window=31,
    figsize=(5, 4),
    standard_scale='var', cmap='RdBu_r',
    use_fitted=True, show_row_names=True, border=False, show=False,
)
```

## Validation

- After `compute_macrostates`: `g.macrostates.cat.categories` should include the cell-type clusters that look like canonical trunks AND terminals. If only one macrostate emerges, `n_states` was too small or the velocity is too noisy.
- After `predict_terminal_states`: terminal categories should *not* include obvious trunks (e.g. 'Ductal' in pancreas). If a trunk shows up as a terminal, the kernel mix is wrong (try more velocity weight) or `n_states` is wrong.
- Fate probabilities in `g.fate_probabilities`: each row sums to ~1 across the `n_terminals` columns. If they don't, the GPCCA didn't converge — try `n_states+1` for macrostates and re-fit.
- After `branch_streamplot`: visually confirm the trunk-to-branch transition aligns with the `branch_center` value. Adjust `branch_center` to position the split point in the figure.
- After `dynamic_features(layer='Ms')`: marker genes peak at biologically expected pseudotimes (e.g. Ins1/Ins2 should peak late on Beta lineage). If they peak early, the pseudotime is reversed — check `velocity_pseudotime` polarity.
- Trunk filtering for `split_time`: the median pseudotime of trunk clusters should be *between* the pre-branch and post-branch peaks. Outside that window, `shared_trunk=True` will look weird.

## Resource Map

- See [`reference.md`](reference.md) for compact copy-paste snippets.
- See [`references/source-grounding.md`](references/source-grounding.md) for which APIs are external (CellRank) vs OmicVerse, and verified signatures for the OmicVerse-side helpers.
- For RNA velocity computation, see the existing `single-cell-rna-velocity` skill.
- For other trajectory backends, see the per-method skills.

## Examples
- "Compute a CellRank fate map on a scvelo-processed pancreas AnnData with `0.8*vk + 0.2*ck` and predict 4 terminal states."
- "Pull the Beta fate probability into `obs['beta_fate']` and plot it on the UMAP alongside cluster labels."
- "Build a branch-aware gene-trend plot for Alpha vs Beta with a shared endocrine trunk."
- "Make a multi-module dynamic heatmap of progenitor + 3 endocrine fates with `use_fitted=True`."

## References
- Tutorial notebook: [`t_velo_cellrank.ipynb`](https://omicverse.readthedocs.io/en/latest/Tutorials-velo/t_velo_cellrank/) — pancreas endocrinogenesis end-to-end with scvelo + CellRank.
- CellRank 2 paper: Lange *et al.* 2022, *Nature Methods* — "CellRank for directed single-cell fate mapping".
- Live API verified — see [`references/source-grounding.md`](references/source-grounding.md).
