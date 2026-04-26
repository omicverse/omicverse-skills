# Paired microbiome × metabolomics — quick commands

```python
import omicverse as ov
import numpy as np
import matplotlib.pyplot as plt
ov.plot_set()

# 1) Real paired data — Franzosa 2019 IBD
adata_mb, adata_mt = ov.micro.fetch_franzosa_ibd_2019(
    data_dir='/scratch/.../franzosa_2019',
)
print('microbes:', adata_mb.shape, 'metabolites:', adata_mt.shape)

# 2) Filter to a tractable size for MMvec
ab_rank = np.argsort(-np.asarray(adata_mb.X).sum(axis=0))[:150]
adata_mb_f = adata_mb[:, ab_rank].copy()
adata_mt_named = adata_mt[:, adata_mt.var['name'].notna()].copy()
var_rank = np.argsort(-np.asarray(adata_mt_named.X).var(axis=0))[:200]
adata_mt_f = adata_mt_named[:, var_rank].copy()
adata_mt_f.var_names = adata_mt_f.var['name'].astype(str).values
print('filtered microbes:', adata_mb_f.shape, 'metabolites:', adata_mt_f.shape)
```

## Spearman triage

```python
spear = ov.micro.paired_spearman(
    adata_mb_f, adata_mt_f,
    clr_microbe=True, log1p_metabolite=True,
    min_prevalence=0.1,
)
print(f'{(spear["fdr_bh"] < 0.05).sum()} pairs at FDR 0.05 / {len(spear)} tested')
print(spear.head(10))   # microbe, metabolite, r, pvalue, fdr_bh
```

## CCA shared modes

```python
cca = ov.micro.paired_cca(
    adata_mb_f, adata_mt_f,
    n_components=3,
    clr_microbe=True, log1p_metabolite=True,
)
print('canonical r:', [round(c, 3) for c in cca['canonical_correlations']])

ov.micro.plot_cca_scatter(cca, component=1)
plt.tight_layout(); plt.show()
```

## MMvec — joint low-rank model

```python
mmvec = ov.micro.MMvec(
    n_latent=4,
    lr=0.05,
    epochs=600,
    val_frac=0.15,
    patience=80,
    l2=1e-3,
    seed=0,
).fit(adata_mb_f, adata_mt_f)
print('best epoch:', mmvec.best_epoch_,
      ' final train loss:', round(mmvec.loss_history_[-1], 4))

# diagnostics
ov.micro.plot_mmvec_training(mmvec)
plt.tight_layout(); plt.show()

# microbe × metabolite log-odds heatmap (top 15 by absolute marginal score)
ov.micro.plot_cooccurrence(mmvec.cooccurrence(), top_n=15)
plt.tight_layout(); plt.show()

# joint embedding biplot
ov.micro.plot_embedding_biplot(mmvec, components=(0, 1), label_top=8)
plt.tight_layout(); plt.show()
```

## MMvec downstream queries

```python
# Top-N pairs by absolute log-odds score
print(mmvec.top_pairs(n=20))               # microbe, metabolite, score

# Symmetric log-odds matrix
M = mmvec.cooccurrence()                   # microbes × metabolites

# Per-microbe predicted metabolite distribution
P = mmvec.conditional_probabilities()       # rows sum to 1
print(P.loc['<microbe_name>'].sort_values(ascending=False).head(10))

# Latent embeddings as DataFrames
print(mmvec.microbe_embeddings_.head())     # rows = microbes, cols = K1..Kn_latent
print(mmvec.metabolite_embeddings_.head())  # rows = metabolites
```

## Synthetic recovery validation

```python
ad_mb_sim, ad_mt_sim, truth = ov.micro.simulate_paired(
    n_samples=30, n_microbes=40, n_metabolites=20,
    n_pairs=5, effect_range=(1.0, 2.0),
    seed=0,
)
res_sp = ov.micro.paired_spearman(ad_mb_sim, ad_mt_sim)
mmvec_sim = ov.micro.MMvec(n_latent=3, epochs=400, val_frac=0.1, seed=0)\
    .fit(ad_mb_sim, ad_mt_sim)

ov.micro.plot_paired_method_comparison(
    truth,
    spearman_df=res_sp,
    mmvec_model=mmvec_sim,
)
plt.tight_layout(); plt.show()
```
