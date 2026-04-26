# Phylogenetic diversity — quick commands

```python
import omicverse as ov
import anndata as ad
import matplotlib.pyplot as plt, pandas as pd
from pathlib import Path
ov.plot_set()

BASE = Path('/scratch/.../cache/16s/run_mothur_sop')
adata = ad.read_h5ad(BASE / 'mothur_sop_16s.h5ad')

# 1) Build the tree from the ASV centroid FASTA
tree = ov.alignment.build_phylogeny(
    asvs_fasta=str(BASE / 'asv' / 'asvs.fasta'),
    workdir=str(BASE / 'phylogeny'),
    mafft_mode='auto',           # 'auto' | 'L-INS-i' | 'E-INS-i'
    fasttree_model='gtr',        # 'gtr' default; 'jc' is didactic only
    gamma=True,
    mafft_threads=8,
    fasttree_threads=4,
)
print('tree path:', tree['tree'])

# 2) Attach to AnnData (prunes tree tips not in adata.var_names by default)
ov.micro.attach_tree(adata, newick=tree['newick'])
print('attached tips:', adata.uns['micro']['tree_tips'])
```

## Faith PD alongside Shannon / Observed

```python
ov.micro.Alpha(adata).run(metrics=['shannon', 'observed_otus', 'faith_pd'])
print(adata.obs[['shannon', 'observed_otus', 'faith_pd']].describe())
```

## UniFrac (unweighted + weighted) + Bray-Curtis baseline

```python
b = ov.micro.Beta(adata)
b.run(metric='unweighted_unifrac', rarefy=False)
b.run(metric='weighted_unifrac',   rarefy=False)
b.run(metric='braycurtis',         rarefy=True)

# Inspect the off-diagonal magnitudes
import numpy as np
for key in ['unweighted_unifrac', 'weighted_unifrac', 'braycurtis']:
    M = adata.obsp[key]
    triu = np.triu_indices_from(M, k=1)
    print(f'{key:20s} mean={M[triu].mean():.3f}  max={M[triu].max():.3f}')
```

## Three-panel PCoA comparison

```python
palette = {'Early': '#2a9d8f', 'Late': '#e76f51', 'Mock': '#264653'}
groups = ['Early', 'Late', 'Mock']

fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
for ax, dist_key, title in zip(
    axes,
    ['braycurtis', 'unweighted_unifrac', 'weighted_unifrac'],
    ['Bray-Curtis', 'unweighted UniFrac', 'weighted UniFrac'],
):
    ord_ = ov.micro.Ordinate(adata, dist_key=dist_key); ord_.pcoa(n=2)
    pct = ord_.proportion_explained() * 100.0
    coords = pd.DataFrame(adata.obsm[f'{dist_key}_pcoa'],
                          index=adata.obs_names, columns=['PC1', 'PC2'])
    for g in groups:
        sub = adata.obs[adata.obs.group == g]
        if sub.empty: continue
        ax.scatter(coords.loc[sub.index, 'PC1'],
                   coords.loc[sub.index, 'PC2'],
                   color=palette[g], s=70, alpha=0.85,
                   edgecolor='k', label=g)
    ax.set_xlabel(f'PC1 ({pct[0]:.1f}%)')
    ax.set_ylabel(f'PC2 ({pct[1]:.1f}%)')
    ax.set_title(title)
axes[-1].legend(title='group', frameon=True)
plt.tight_layout(); plt.show()
```

## Strict tip-name validation

```python
# Raises if any tree tip is not in adata.var_names (use in CI)
ov.micro.attach_tree(adata, newick=tree['newick'], strict=True)
```

## Persist

```python
adata.write_h5ad(BASE / 'mothur_sop_16s_tree.h5ad')
```
