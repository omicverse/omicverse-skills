# Bulk cell-type deconvolution — quick commands

## Setup — load bulk + single-cell reference

```python
import omicverse as ov
ov.plot_set()

# Built-in COVID-19 PBMC demo (replace with your own)
bulk_ad        = ov.datasets.decov_bulk_covid_bulk()
single_ad_ref  = ov.datasets.decov_bulk_covid_single()

# Reference QC
ov.pl.embedding(single_ad_ref, basis='X_umap',
                color=['cell.type.fine', 'cell.type.coarse'], ncols=1)
ov.pl.embedding(single_ad_ref, basis='X_umap', color='disease')

fig, ax = ov.plt.subplots(figsize=(1.5, 3))
ov.pl.cellproportion(
    adata=single_ad_ref,
    celltype_clusters='cell.type.coarse',
    groupby='disease', legend=True, ax=ax,
)
```

## BayesPrism (multi-core CPU, full Bayesian)

```python
deconv_obj = ov.bulk.Deconvolution(
    adata_bulk=bulk_ad,
    adata_single=single_ad_ref,
    max_single_cells=10000,
    celltype_key='cell.type.coarse',
    cellstate_key='cell.type.fine',          # hierarchical signature
)
res = deconv_obj.deconvolution(
    method='bayesprism',
    n_cores=8, fast_mode=True,
)
res = res[single_ad_ref.obs['cell.type.coarse'].cat.categories]
print(res.head())
```

## Scaden (DNN, GPU-friendly)

```python
deconv_obj = ov.bulk.Deconvolution(
    adata_bulk=bulk_ad,
    adata_single=single_ad_ref,
    max_single_cells=10000,
    celltype_key='cell.type.coarse',
    cellstate_key='cell.type.fine',
    gpu='mps',                               # or gpu=0 for CUDA, gpu=-1 for CPU
)
res2 = deconv_obj.deconvolution(
    method='scaden',
    scaler='ss', scale=True, datatype='counts',
    pseudobulk_size=2000,
)
res2 = res2[single_ad_ref.obs['cell.type.coarse'].cat.categories]
```

## TAPE (default; CPU-friendly autoencoder)

```python
deconv_obj = ov.bulk.Deconvolution(
    adata_bulk=bulk_ad,
    adata_single=single_ad_ref,
    max_single_cells=10000,
    celltype_key='cell.type.coarse',
)
res_tape = deconv_obj.deconvolution(
    method='tape',
    scaler='mms',                            # TAPE default
    datatype='counts',
    adaptive=True, variance_threshold=0.98,
    batch_size=128, epochs=128, seed=1,
)
```

## OmicsTweezer

```python
res_ot = deconv_obj.deconvolution(method='omicstweezer')
```

## Stacked-bar visualisation

```python
import matplotlib.pyplot as plt

color_dict = dict(zip(
    single_ad_ref.obs['cell.type.coarse'].cat.categories,
    single_ad_ref.uns['cell.type.coarse_colors'],
))

ax = res.plot(kind='bar', stacked=True, figsize=(12, 3), color=color_dict)
ax.set(xlabel='Sample', ylabel='Cell Fraction',
       title='BayesPrism predicted fractions')
ov.plt.legend(bbox_to_anchor=(1.05, 1), ncol=1)
ov.plt.show()
```

## Group-averaged stacked bars

```python
ov.pl.plot_grouped_fractions(
    res,
    obs=bulk_ad.obs,
    group_key='disease',                     # phenotype column on bulk side
    color_dict=color_dict,
    agg='mean',
    normalize=True,
    figsize=(4, 4),
)
```

## Cross-method agreement check

```python
import numpy as np

common = res.columns.intersection(res2.columns)
per_sample_corr = np.array([
    np.corrcoef(res.loc[s, common].values,
                res2.loc[s, common].values)[0, 1]
    for s in res.index
])
print(f'BayesPrism vs Scaden per-sample r: '
      f'{per_sample_corr.mean():.3f} ± {per_sample_corr.std():.3f}')
print(f'Samples with r < 0.5: '
      f'{(per_sample_corr < 0.5).sum()} / {len(per_sample_corr)}')
```

## Method matrix (when to pick which)

| method        | Engine    | Compute       | Use when                                     |
|---------------|-----------|---------------|----------------------------------------------|
| `tape`        | DNN (AE)  | CPU OK        | Default; sensible on most cohorts            |
| `scaden`      | DNN       | GPU helps     | Diverse cohorts; per-cohort training         |
| `bayesprism`  | Bayesian  | Multi-core CPU | Posterior fractions; small / noisy cohorts; needs `cellstate_key` for hierarchy |
| `omicstweezer`| Joint     | CPU OK        | Reference-correction integrated              |
