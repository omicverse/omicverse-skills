# MetaTiME annotation — quick commands

```python
import omicverse as ov
import scanpy as sc
ov.utils.ov_plot_set()

# 1) Load batch-corrected AnnData
adata = sc.read('TiME_adata_scvi.h5ad')

# 2) MDE projection (optional, faster than UMAP)
adata.obsm['X_mde'] = ov.utils.mde(adata.obsm['X_scVI'])

# 3) Three-step MetaTiME
TiME_object = ov.single.MetaTiME(adata, mode='table')
TiME_object.overcluster(resolution=8, clustercol='overcluster')
TiME_object.predictTiME(save_obs_name='MetaTiME')

# obs columns added: 'overcluster', 'MetaTiME', 'Major_MetaTiME'
print(adata.obs[['MetaTiME', 'Major_MetaTiME']].head())
print(adata.obs['Major_MetaTiME'].value_counts())
```

## Built-in labelled plot

```python
fig, ax = TiME_object.plot(
    cluster_key='MetaTiME',
    basis='X_mde',
    fontsize=8,
    min_cell=5,                # label only clusters with >= 5 cells
    figsize=(6, 6),
    dpi=80,
)
```

## Standard scanpy embedding for major categories

```python
sc.pl.embedding(
    adata, basis='X_mde',
    color=['Major_MetaTiME'],
    frameon=False, ncols=1,
)

# Or both fine + major side-by-side
sc.pl.embedding(
    adata, basis='X_mde',
    color=['MetaTiME', 'Major_MetaTiME'],
    frameon=False, ncols=2, wspace=0.4,
)
```

## Adjusting over-clustering resolution

```python
# For a >100k-cell cohort, push resolution higher
TiME_big = ov.single.MetaTiME(adata, mode='table')
TiME_big.overcluster(resolution=10)
TiME_big.predictTiME(save_obs_name='MetaTiME_high_res')

# For a <5k-cell cohort, lower resolution
TiME_small = ov.single.MetaTiME(adata_small, mode='table')
TiME_small.overcluster(resolution=4)
TiME_small.predictTiME(save_obs_name='MetaTiME')
```

## Inspect MeC scoring distribution

```python
# adata.obs['MetaTiME'] is the fine label; cell-state distribution
fine_counts = adata.obs.groupby('Major_MetaTiME')['MetaTiME'].value_counts()
print(fine_counts.head(20))
```
