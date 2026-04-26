# CellRank fate maps — quick commands

```python
import omicverse as ov
import scanpy as sc
import scvelo as scv
import cellrank as cr
import numpy as np
import matplotlib.pyplot as plt
ov.plot_set(font_path='Arial')

# Assumes velocity is computed (see single-cell-rna-velocity skill):
# adata.layers['velocity'], uns['velocity_graph'],
# obs['development_pseudotime'], layers['Ms'], obs['clusters']

# 1) Kernel mix: velocity drives direction, connectivity smooths noise
vk = cr.kernels.VelocityKernel(adata).compute_transition_matrix()
ck = cr.kernels.ConnectivityKernel(adata).compute_transition_matrix()
kernel = 0.8 * vk + 0.2 * ck

# 2) GPCCA → macrostates → terminal states
g = cr.estimators.GPCCA(kernel)
g.compute_macrostates(n_states=6, cluster_key='clusters')
print('macrostates:', list(g.macrostates.cat.categories))

g.predict_terminal_states(n_states=4)
print('terminals:', list(g.terminal_states.cat.categories))

# 3) Per-cell fate probability into adata.obs
adata.obs['beta_fate'] = np.asarray(g.fate_probabilities['Beta']).ravel()

ov.pl.embedding(adata, basis='X_umap',
                color=['clusters', 'beta_fate'],
                cmap='Reds', frameon='small')
```

## Branch streamplot (pseudotime + clusters)

```python
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

## Per-fate gene trends (Beta lineage)

```python
adata_beta = adata[adata.obs['beta_fate'] > 0.15].copy()

beta_dyn = ov.single.dynamic_features(
    adata_beta,
    genes=['Pdx1', 'Ins1', 'Ins2', 'Sox9'],
    pseudotime='development_pseudotime',
    layer='Ms',                 # smoothed-spliced from scvelo
    distribution='normal', link='identity',
    n_splines=8,
    store_raw=True, raw_obs_keys=['clusters'],
)
ov.pl.dynamic_trends(
    beta_dyn, genes=['Pdx1', 'Ins1', 'Ins2', 'Sox9'],
    compare_features=True,
    add_point=True, point_color_by='clusters',
    line_style_by='features',
    figsize=(6.2, 3.2), linewidth=2.2,
    legend_loc='right margin', legend_fontsize=8,
)
```

## Branch-aware comparison (Alpha vs Beta with shared endocrine trunk)

```python
branch_clusters = ['Alpha', 'Beta']
split_mask = adata.obs['clusters'].astype(str).isin(
    ['Ngn3 high EP', 'Pre-endocrine']
)

cr_branch_dyn = ov.single.dynamic_features(
    adata,
    genes=['Gcg', 'Ins1', 'Ins2', 'Pdx1'],
    pseudotime='development_pseudotime',
    layer='Ms',
    groupby='clusters', groups=branch_clusters,
    distribution='normal', link='identity',
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

## Multi-module dynamic heatmap

```python
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
