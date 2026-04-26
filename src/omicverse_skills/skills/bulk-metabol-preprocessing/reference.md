# Metabolomics preprocessing — quick commands

## Functional API (one-shot)

```python
import omicverse as ov
import numpy as np

ov.plot_set()

# 1) load — replace with your file
csv_path = ov.datasets.download_data(
    url='https://rest.xialab.ca/api/download/metaboanalyst/human_cachexia.csv',
    file_path='human_cachexia.csv',
    dir='./metabol_data',
)
adata = ov.metabol.read_metaboanalyst(csv_path, group_col='Muscle loss')

# 2) imputation (qrilc for LC-MS / left-censored; knn for MAR)
adata = ov.metabol.impute(adata, method='qrilc', q=0.01, seed=0)

# 3) normalize per sample (PQN canonical)
adata = ov.metabol.normalize(adata, method='pqn')

# 4) transform per feature
log_adata = ov.metabol.transform(adata, method='log')               # for univariate stats
pareto_adata = ov.metabol.transform(log_adata, method='pareto',     # for PLS-DA / OPLS-DA
                                    stash_raw=False)

# 5) downstream stats use log_adata
deg = ov.metabol.differential(
    log_adata,
    group_col='group', group_a='cachexic', group_b='control',
    method='welch_t', log_transformed=True,
)
fig, ax = ov.metabol.volcano(deg, padj_thresh=0.10, log2fc_thresh=0.3, label_top_n=8)
```

## Chained class API (`pyMetabo`)

```python
m = (
    ov.metabol.pyMetabo(adata.copy())
      .impute(method='qrilc', seed=0)
      .normalize(method='pqn')
      .transform(method='log')
      .differential(group_col='group', group_a='cachexic', group_b='control',
                    method='welch_t', log_transformed=True)
      .transform(method='pareto', stash_raw=False)
      .opls_da(n_ortho=1)
)
m.deg_table.head()
m.vip_table().head()
m.significant_metabolites(padj_thresh=0.10, log2fc_thresh=0.3)
```

## LC-MS drift / SERRF / ComBat

```python
# adata.obs needs columns: 'sample_type' (QC|real), 'batch', 'order'

# 1) injection-order drift (LOESS on QC samples)
adata = ov.metabol.drift_correct(
    adata,
    injection_order='order',
    qc_mask=(adata.obs['sample_type'] == 'QC').to_numpy(),
    frac=0.5,
)

# 2) SERRF — Random Forest residual correction trained on QC pools
adata = ov.metabol.serrf(
    adata,
    qc_col='sample_type', qc_label='QC',
    batch_col='batch',
    top_k=8, n_estimators=100, seed=0,
)
# inspect: adata.var[['cv_qc_raw', 'cv_qc_serrf']]

# 3) ComBat — between-batch additive shifts (run AFTER serrf, only if residual batch effect)
ov.bulk.batch_correction(adata, batch_key='batch', key_added='batch_correction')
# corrected matrix at adata.layers['batch_correction']
```

## Sample-level outlier QC

```python
real = adata[adata.obs['sample_type'] == 'real'].copy()
qc_df = ov.metabol.sample_qc(real, n_components=3, alpha=0.95)
ov.metabol.sample_qc_plot(qc_df)
flagged = qc_df.index[qc_df['is_outlier']].tolist()
```

## Multi-group ANOVA (3+ groups)

```python
aov = ov.metabol.anova(adata, group_col='dose_tertile', method='welch_anova')
# also supports method='anova' (assumes equal var) or 'kruskal' (non-parametric)
print((aov['padj'] < 0.05).sum(), 'significant features')
```
