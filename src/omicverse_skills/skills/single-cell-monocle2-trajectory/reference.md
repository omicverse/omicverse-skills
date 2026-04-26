# Monocle2 trajectory — quick commands

```python
import anndata as ad
import omicverse as ov
from omicverse.single import Monocle
import matplotlib.pyplot as plt
ov.plot_set(font_path='Arial')

# 1) Load preprocessed AnnData (e.g. Olsson WT subset)
adata = ad.read_h5ad('olsson_wt.h5ad')

# 2) Fit
mono = Monocle(adata)
mono.preprocess()
mono.select_ordering_genes(max_genes=1000)
mono.plot_ordering_genes(figsize=(5, 4)); plt.show()

# 3) DDRTree trajectory (auto-fits on first plot call)
mono.plot_trajectory(
    color_by='State', cell_size=2.4, cell_link_size=0.5,
    show_branch_points=False, figsize=(6, 4),
); plt.show()

# 4) Branch streamplot view of pseudotime by celltype
fig, ax = ov.pl.branch_streamplot(
    mono.adata, group_key='subtype',
    pseudotime_key='Pseudotime', show=False,
); plt.show()
```

## DEG along pseudotime + heatmap

```python
ordering = mono.adata.var_names[mono.adata.var['use_for_ordering']].tolist()
mono_ord = Monocle(mono.adata[:, ordering].copy())   # subset for speed

de = mono_ord.differential_gene_test(cores=-1)
sig = de[(de['qval'] < 0.01) & (de['status'] == 'OK')]
top40 = sig.sort_values('pval').head(40).index.tolist()

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
```

## Per-gene GAM trends

```python
markers = [g for g in ['Gfi1', 'Irf8', 'Elane', 'Prtn3', 'Mpo', 'Car2']
           if g in mono.adata.var_names]

res = ov.single.dynamic_features(
    mono.adata,
    genes=markers,
    pseudotime='Pseudotime',
    distribution='normal', link='identity',
    n_splines=8, spline_order=3,
    store_raw=True, raw_obs_keys=['State'],
)

ov.pl.dynamic_trends(
    res, genes=markers,
    add_point=True, point_color_by='State',
    figsize=(4, 3.5),
    legend_loc='right margin', legend_fontsize=8,
); plt.show()
```

## Compare features on the same axis

```python
trend_genes = ['Gfi1', 'Irf8', 'Elane', 'Car2']
res_global = ov.single.dynamic_features(
    mono.adata, genes=trend_genes, pseudotime='Pseudotime',
    store_raw=True, raw_obs_keys=['subtype', 'State'],
)
ov.pl.dynamic_trends(
    res_global, genes=trend_genes,
    compare_features=True,
    add_point=True, point_color_by='subtype',
    line_style_by='features',
    figsize=(6, 3.2), linewidth=2.2,
)
```

## Branch-aware trends with a shared trunk

```python
import numpy as np

branch_subtypes = ['Gmp', 'LK']     # adjust to your trajectory's branches
res_branch = ov.single.dynamic_features(
    mono.adata,
    genes=trend_genes,
    pseudotime='Pseudotime',
    groupby='subtype',
    groups=branch_subtypes,
    store_raw=True,
)

# Pick the median pseudotime of the trunk subtype as the split time
trunk_mask = mono.adata.obs['subtype'].astype(str).isin(['Cmp'])
split_time = float(np.nanmedian(mono.adata.obs.loc[trunk_mask, 'Pseudotime']))

ov.pl.dynamic_trends(
    res_branch, genes=trend_genes,
    compare_groups=True,
    split_time=split_time,
    shared_trunk=True,
    add_point=True, point_color_by='group',
    figsize=(4.6, 3), ncols=2,
    title='Branch-aware marker trends',
)
```

## BEAM branch-dependent genes

```python
beam = mono_ord.BEAM(branch_point=1, cores=-1)
sig_beam = beam[beam['qval'] < 0.01].sort_values('qval')
print(f'BEAM hits: {len(sig_beam)} / {len(beam)}')

top_branch = sig_beam.head(4).index.tolist()
beam_dyn = ov.single.dynamic_features(
    mono.adata,
    genes=top_branch,
    pseudotime='Pseudotime',
    groupby='subtype', groups=['Gmp', 'LK'],
    store_raw=True,
)
ov.pl.dynamic_trends(
    beam_dyn, genes=top_branch,
    compare_groups=True,
    split_time=split_time, shared_trunk=True,
    add_point=True, point_color_by='group',
    figsize=(6, 4), ncols=2, linewidth=2.2,
    title='BEAM branch-dependent genes',
)
```

## Multi-set heatmap (gene panels)

```python
panels = {
    'Early progenitor': ['Gfi1'],
    'Granulocyte-like': ['Elane', 'Prtn3', 'Mpo'],
    'Alternative':       ['Irf8', 'Car2'],
}
panels = {k: [g for g in v if g in mono.adata.var_names] for k, v in panels.items()}

g = ov.pl.dynamic_heatmap(
    mono.adata,
    var_names=panels,
    pseudotime='Pseudotime',
    cell_annotation='subtype',
    cell_bins=140, smooth_window=13, fitted_window=25,
    figsize=(5, 4),
    standard_scale='var', cmap='RdBu_r',
    use_fitted=True, border=False, show=False,
)
```
