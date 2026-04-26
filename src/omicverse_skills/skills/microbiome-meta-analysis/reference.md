# Cross-cohort meta-analysis — quick commands

```python
import omicverse as ov
import anndata as ad
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
ov.plot_set()

# Per-cohort AnnDatas (each from the 16S amplicon skill)
adata_A = ad.read_h5ad('cohort_A_genus.h5ad')
adata_B = ad.read_h5ad('cohort_B_genus.h5ad')
adata_C = ad.read_h5ad('cohort_C_genus.h5ad')

# 1) Combine — collapse to genus, zero-fill cohort-specific features, tag obs['study']
adata_all = ov.micro.combine_studies(
    [adata_A, adata_B, adata_C],
    study_names=['cohort_A', 'cohort_B', 'cohort_C'],
    rank='genus',
)
print('combined:', adata_all.shape)
print(adata_all.obs['study'].value_counts())

# 2) Diagnostic PCoA on the combined object
ov.micro.Beta(adata_all).run(metric='braycurtis', rarefy=False)
coords = np.asarray(ov.micro.Ordinate(adata_all, dist_key='braycurtis').pcoa(n=2))

palette = {'cohort_A': '#E57373', 'cohort_B': '#64B5F6', 'cohort_C': '#81C784'}
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
for st, color in palette.items():
    m = adata_all.obs['study'].values == st
    axes[0].scatter(coords[m, 0], coords[m, 1], c=color, label=st,
                    s=40, edgecolors='k', linewidths=0.5, alpha=0.85)
axes[0].set(xlabel='PCo1', ylabel='PCo2', title='Coloured by cohort')
axes[0].legend()

for grp, marker, c in [('CTRL', 'o', '#BDBDBD'), ('CASE', '^', '#7E57C2')]:
    m = adata_all.obs['disease'].values == grp
    axes[1].scatter(coords[m, 0], coords[m, 1], marker=marker, c=c, label=grp,
                    s=40, edgecolors='k', linewidths=0.5, alpha=0.85)
axes[1].set(xlabel='PCo1', ylabel='PCo2', title='Coloured by disease')
axes[1].legend()
plt.tight_layout(); plt.show()
```

## Per-study DA (sanity check)

```python
per_study_top = {}
for name, a in [('A', adata_A), ('B', adata_B), ('C', adata_C)]:
    a_g = ov.micro.collapse_taxa(a, rank='genus')
    res = ov.micro.DA(a_g).deseq2(
        group_key='disease', group_a='CTRL', group_b='CASE',
        min_prevalence=0.1,
    )
    per_study_top[name] = res.head(5)['feature'].tolist()
    print(f'{name} top-5: {per_study_top[name]}')
print('overlap:', set.intersection(*map(set, per_study_top.values())))
```

## Random-effects meta-analysis

```python
meta = ov.micro.meta_da(
    [adata_A, adata_B, adata_C],
    study_names=['cohort_A', 'cohort_B', 'cohort_C'],
    group_key='disease', group_a='CTRL', group_b='CASE',
    method='deseq2',                  # also: 'wilcoxon', 'ancombc'
    rank='genus',
    min_prevalence=0.1,
    combine='random_effects',         # 'fixed_effects' only when I² is provably zero
)

cols = ['feature', 'combined_lfc', 'combined_se', 'z',
        'p_value', 'fdr_bh', 'n_studies', 'I2']
print(f'tested cross-cohort: {len(meta)};  significant at FDR 0.05: {int((meta["fdr_bh"]<0.05).sum())}')
print(meta[cols].head(10))
```

## Forest-style plot of top features

```python
top = meta.head(12).copy().iloc[::-1]
ci = 1.96 * top['combined_se'].values
fig, ax = plt.subplots(figsize=(7, 5))
ax.barh(top['feature'], top['combined_lfc'], xerr=ci,
        edgecolor='black', linewidth=0.4)
ax.axvline(0, color='k', lw=0.7)
ax.set_xlabel('combined log2 fold-change (CASE / CTRL)')
plt.tight_layout(); plt.show()
```

## Effect × heterogeneity diagnostic

```python
sig = meta[meta['fdr_bh'] < 0.05].copy()
fig, ax = plt.subplots(figsize=(7, 4))
ax.scatter(sig['combined_lfc'], sig['I2'], s=50, edgecolors='k', linewidths=0.4)
ax.axhline(0.25, color='grey', lw=0.5, ls='--')   # low / moderate I²
ax.axhline(0.75, color='grey', lw=0.5, ls='--')   # substantial / high I²
ax.set_xlabel('combined log2FC'); ax.set_ylabel("Cochran's I²")
plt.tight_layout(); plt.show()
```

## Method-cycling benchmark

```python
results = {
    method: ov.micro.meta_da(
        [adata_A, adata_B, adata_C],
        study_names=['cohort_A', 'cohort_B', 'cohort_C'],
        group_key='disease', group_a='CTRL', group_b='CASE',
        method=method, rank='genus', min_prevalence=0.1,
    )
    for method in ['wilcoxon', 'deseq2', 'ancombc']
}

summary = pd.DataFrame({
    'n_tested':       {m: len(r) for m, r in results.items()},
    'n_sig_fdr_0.05': {m: int((r['fdr_bh'] < 0.05).sum()) for m, r in results.items()},
})
print(summary)
```
