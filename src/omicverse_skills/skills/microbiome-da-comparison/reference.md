# DA method comparison — quick commands

```python
import omicverse as ov
import anndata as ad
import matplotlib.pyplot as plt

ov.plot_set()

# Reuse AnnData from the 16S amplicon skill; restrict to a clean two-group contrast
adata = ad.read_h5ad('mothur_sop_16s.h5ad')
adata = adata[adata.obs['group'].isin(['Early', 'Late'])].copy()

# Collapse to genus before DA — the canonical reporting rank for 16S
adata_genus = ov.micro.collapse_taxa(adata, rank='genus')

# Run all three on the same contrast
wx = ov.micro.DA(adata_genus).wilcoxon(
    group_key='group', group_a='Early', group_b='Late', min_prevalence=0.1,
)
ds = ov.micro.DA(adata_genus).deseq2(
    group_key='group', group_a='Early', group_b='Late', min_prevalence=0.1,
)
ab = ov.micro.DA(adata_genus).ancombc(
    group_key='group', min_prevalence=0.1,
)

# Column-name pitfall: ANCOM-BC may report q_value or fdr_bh depending on skbio version
sig_col_ab = 'q_value' if 'q_value' in ab.columns else 'fdr_bh'

# Build hit sets at FDR 0.05
sig_wx = set(wx.loc[wx['fdr_bh']  < 0.05, 'feature'])
sig_ds = set(ds.loc[ds['fdr_bh']  < 0.05, 'feature'])
sig_ab = set(ab.loc[ab[sig_col_ab] < 0.05, 'feature'])

print('Wilcoxon only :', len(sig_wx - sig_ds - sig_ab))
print('DESeq2 only   :', len(sig_ds - sig_wx - sig_ab))
print('ANCOM-BC only :', len(sig_ab - sig_wx - sig_ds))
print('Wx ∩ DESeq2   :', len((sig_wx & sig_ds) - sig_ab))
print('Wx ∩ ANCOM-BC :', len((sig_wx & sig_ab) - sig_ds))
print('DS ∩ ANCOM-BC :', len((sig_ds & sig_ab) - sig_wx))
print('All three     :', len(sig_wx & sig_ds & sig_ab))

# Optional Venn diagram
try:
    from matplotlib_venn import venn3
    fig, ax = plt.subplots(figsize=(5, 5))
    venn3([sig_wx, sig_ds, sig_ab],
          set_labels=('Wilcoxon', 'DESeq2', 'ANCOM-BC'), ax=ax)
    ax.set_title('Genera significant at FDR/q < 0.05')
    plt.tight_layout(); plt.show()
except ImportError:
    print('(install matplotlib_venn for a rendered Venn)')
```

## Reporting strategies

```python
consensus    = sig_wx & sig_ds & sig_ab               # strongest claim
defensible   = sig_wx & sig_ab                         # bypasses NB assumption, keeps compositional + non-parametric
exploratory  = sig_wx | sig_ds | sig_ab                # union — flag method-specific hits explicitly

print('CONSENSUS:', sorted(consensus))
print('DEFENSIBLE (Wx ∩ ANCOM-BC):', sorted(defensible))
print('EXPLORATORY (any method):', len(exploratory))
```

## Disagreement diagnostic

```python
# Find DESeq2-only hits — likely NB-shrinkage artefacts on rare features
deseq2_only = sig_ds - sig_wx - sig_ab
if deseq2_only:
    suspect = adata_genus.var.copy()
    counts = adata_genus.X.toarray() if hasattr(adata_genus.X, 'toarray') else adata_genus.X
    suspect['n_nonzero'] = (counts > 0).sum(axis=0)
    suspect['mean_count'] = counts.mean(axis=0)
    print(suspect.loc[list(deseq2_only), ['n_nonzero', 'mean_count']])
```
