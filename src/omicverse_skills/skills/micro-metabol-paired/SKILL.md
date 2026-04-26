---
name: omicverse-micro-metabol-paired
description: Paired microbiome × metabolomics integration on sample-aligned AnnDatas. Provides three methods - Spearman correlation (fast pairwise FDR), sklearn CCA (linear canonical mode shared between modalities), and MMvec (PyTorch low-rank co-occurrence model with conditional probabilities and biplot embeddings). Use when you have paired microbe-counts + metabolite-intensity tables on the same samples and want to find microbe-metabolite pairs that covary, or build a shared latent space.
---

# OmicVerse Microbiome × Metabolomics — Paired Integration

## Goal

Take two `AnnData` objects sharing the same `obs_names` in the same order — one with microbe / ASV / OTU counts, one with metabolite intensities — and find microbe-metabolite pairs that co-vary across samples. Three methods, three different statistical assumptions:

1. **Spearman pairwise correlation** (`paired_spearman`) — fast; per-pair `r` + p + BH-FDR; symmetric / non-parametric.
2. **CCA (`paired_cca`)** — linear canonical components; finds **shared latent modes**; small number of components, each explaining a covariance pair.
3. **MMvec (`MMvec`)** — PyTorch low-rank model `logits = U·V.T + β` over per-microbe metabolite distributions; gives an **embedding biplot** + **conditional probabilities** P(metabolite | microbe).

Plus a synthetic-recovery validation harness (`simulate_paired`) and a method-comparison plot (`plot_paired_method_comparison`).

This is a **two-modality integration** skill; the upstream microbe table is presumed to come from the 16S amplicon skill, and the metabolite table from any of the metabolomics skills. The Franzosa 2019 IBD dataset is the bundled real-data demo (`fetch_franzosa_ibd_2019`).

## Quick Workflow

1. Load or build the two paired AnnDatas. They must share `obs_names` in the same order — the methods enforce this.
2. **Filter to a tractable size**: MMvec scales `M × N × K` (microbes × metabolites × latent dim); the Spearman matrix is `M × N` pairs. For Franzosa: top 150 microbes by total counts, top 200 metabolites by variance is the canonical demo trim.
3. **Spearman first** — it's free: `ov.micro.paired_spearman(adata_mb, adata_mt, clr_microbe=True, log1p_metabolite=True, min_prevalence=0.1) → pd.DataFrame` with `microbe`, `metabolite`, `r`, `pvalue`, `fdr_bh`. Use the FDR cutoff to triage candidate pairs.
4. **CCA** for shared linear modes: `cca = ov.micro.paired_cca(adata_mb, adata_mt, n_components=3, clr_microbe=True, log1p_metabolite=True) → dict`. Inspect `cca['canonical_correlations']` — values close to 1 mean strong shared mode; close to 0 means no shared linear signal.
5. **Plot CCA scatter**: `ov.micro.plot_cca_scatter(cca, component=1)` shows the canonical scores from microbiome side vs metabolomics side along the requested component.
6. **MMvec** for richer co-occurrence: `mmvec = ov.micro.MMvec(n_latent=4, epochs=600, val_frac=0.15, patience=80, seed=0).fit(adata_mb, adata_mt)`. `best_epoch_` reports the early-stop checkpoint.
7. **Inspect MMvec output**: `plot_mmvec_training(mmvec)` (loss curve), `plot_cooccurrence(mmvec.cooccurrence(), top_n=15)` (microbe × metabolite log-odds heatmap), `plot_embedding_biplot(mmvec, components=(0, 1), label_top=8)` (joint biplot in latent space). For specific queries use `mmvec.top_pairs(n=20)` and `mmvec.conditional_probabilities()`.
8. **Validate on synthetics**: `adata_sim_mb, adata_sim_mt, truth = simulate_paired(n_pairs=5, seed=0)`; run all three methods on the synthetic data; `plot_paired_method_comparison(truth, spearman_df=..., mmvec_model=...)` shows how close each method got to recovering the planted pairs.

## Interface Summary

Data ingest (real demo + synthetic):
- `ov.micro.fetch_franzosa_ibd_2019(data_dir, overwrite=False, microbe_count_scale=1e6) → (adata_mb, adata_mt)` — auto-downloads + parses the Franzosa 2019 IBD paired dataset (Crohn / UC / Control). Two AnnDatas with identical `obs_names`. ~165 samples × ~600 microbes × ~7000 metabolites by default.
- `ov.micro.simulate_paired(n_samples=30, n_microbes=40, n_metabolites=20, n_pairs=5, effect_range=(1.0, 2.0), depth_range=(1000, 10000), seed=0) → (adata_mb, adata_mt, truth_df)` — synthetic data with `n_pairs` planted producer-pairs.

Pairwise methods:
- `ov.micro.paired_spearman(adata_microbe, adata_metabolite, *, clr_microbe=True, log1p_metabolite=True, min_prevalence=0.0) → pd.DataFrame` columns `microbe`, `metabolite`, `r`, `pvalue`, `fdr_bh`.
- `ov.micro.paired_cca(adata_microbe, adata_metabolite, *, n_components=3, clr_microbe=True, log1p_metabolite=True, max_iter=500) → Dict[str, Any]` with `canonical_correlations`, `microbe_scores`, `metabolite_scores`, `microbe_loadings`, `metabolite_loadings`.

MMvec:
- `MMvec(n_latent=3, lr=0.05, epochs=1000, val_frac=0.1, patience=100, l2=1e-3, seed=0, device=None)` — constructor.
- `.fit(adata_microbe, adata_metabolite, verbose=False) → self` — trains; populates `U_`, `V_`, `beta_`, `microbe_names_`, `metabolite_names_`, `loss_history_`, `val_loss_history_`, `best_epoch_`.
- `.cooccurrence() → pd.DataFrame` — symmetric log-odds matrix `U @ V.T` (microbes × metabolites).
- `.conditional_probabilities() → pd.DataFrame` — softmax-per-row P(metabolite | microbe).
- `.top_pairs(n=20) → pd.DataFrame` columns `microbe`, `metabolite`, `score` — sorted by `|log-odds|`.
- `.microbe_embeddings_` / `.metabolite_embeddings_` — latent-factor DataFrames keyed by name.

Plotting helpers:
- `ov.micro.plot_mmvec_training(mmvec, ax=None)` — loss curve with best-validation epoch marker.
- `ov.micro.plot_cca_scatter(cca_result, component=1, ax=None)` — microbe vs metabolite canonical scores on one component.
- `ov.micro.plot_cooccurrence(score_df, top_n=15, ax=None, cmap='RdBu_r')` — heatmap of a microbe × metabolite score matrix.
- `ov.micro.plot_embedding_biplot(mmvec, components=(0, 1), label_top=5, ax=None)` — joint biplot of microbe + metabolite embeddings in latent space.
- `ov.micro.plot_paired_method_comparison(truth, spearman_df=None, mmvec_model=None, ax=None)` — grouped bar chart of planted-pair recovery rank under each method (works on `simulate_paired` output).

## Boundary

**Inside scope:**
- Sample-aligned microbe + metabolite paired analysis on AnnDatas.
- Three methods (Spearman / CCA / MMvec) on the same input.
- Synthetic-recovery validation harness.
- Franzosa 2019 IBD as the canonical real demo.

**Outside scope — separate skill:**
- Building the per-modality AnnDatas — see `omicverse-microbiome-16s-amplicon-dada2` and `omicverse-bulk-metabol-preprocessing` (or any other metabolomics skill).
- 3+ modality joint factorization (microbe + metabolite + RNA) — that's MOFA territory; see `omicverse-bulk-metabol-pathway-multifactor` for `run_mofa` (which can absorb a microbe view as a third AnnData).
- Causal / mediation analysis between microbes and metabolites — out of scope; `paired_*` reports association only.
- Microbe ↔ host-gene-expression integration — out of scope; same MMvec / CCA shape would work but isn't bundled.

## Branch Selection

**Method choice — by question and cohort size**
- **Need fast triage / FDR-controlled pairs**: Spearman. O(M·N) cost; runs in seconds even on 1000 microbes × 1000 metabolites; gives a per-pair FDR. Best when the question is "give me the top N pairs to chase".
- **Need a small set of interpretable shared modes**: CCA. 2–5 components; each component is a covariance-paired linear axis. Best when reviewers want a "PC1 × PC1" plot and a small canonical-correlation number.
- **Need richer per-microbe distributions / embedding biplot**: MMvec. Slowest (PyTorch training, minutes); gives `P(metabolite | microbe)` and a joint embedding. Best when downstream analysis is "for this microbe, what's the predicted metabolite profile?" or when you want a biplot for a figure.
- **Always run Spearman first** — it doesn't cost anything and gives you a sanity reference for the more elaborate methods.

**`clr_microbe=True` / `log1p_metabolite=True`** — sensible defaults for both Spearman and CCA. CLR on counts handles compositional bias; log1p on metabolites stabilises variance across orders of magnitude. Set to `False` only when the data is already on the appropriate scale (e.g. you've run `ov.metabol.transform(method='log')` already).

**MMvec hyperparameters**
- `n_latent`: 3–5 is canonical. <3 underfits; >10 overfits on small cohorts and slows training quadratically.
- `epochs=600` + `patience=80`: tutorial defaults; the early-stop on validation loss usually triggers well before `epochs`.
- `val_frac=0.15`: 10–20 % validation split. With <30 samples, use `val_frac=0.0` and rely on `epochs` budget alone.
- `lr=0.05` + `l2=1e-3`: stable for the Franzosa-scale demo. Reduce `lr` to 0.01 if loss curves look noisy.
- `seed`: always set for reproducibility.

**Filtering strategy before MMvec**
- Microbes: top 100–200 by total reads (or prevalence). MMvec memory and compute scale linearly with `M` per epoch.
- Metabolites: keep only annotated (filter out unmapped clusters), then top 100–500 by variance. Variance filter is critical because constant or near-constant metabolites contribute no signal but blow training time.

**Reading the cooccurrence heatmap**
- Red cells (positive log-odds) — microbe and metabolite covary positively across samples.
- Blue cells — anti-covariation (rare in real data; usually a normalisation artefact).
- The matrix is dense by construction; `plot_cooccurrence(top_n=15)` truncates to the top-N microbes and top-N metabolites by absolute marginal score so the plot stays legible.

**Conditional probabilities** are useful when the question is asymmetric ("for this microbe, predict the metabolite profile"). The cooccurrence is symmetric; the conditional probabilities are not (each row sums to 1).

## Input Contract

- Both AnnDatas: identical `obs_names` in the same order (`_check_paired` enforces). Length = number of paired samples.
- `adata_microbe.X`: integer counts (or float counts surviving the CLR; the wrappers add a pseudocount internally).
- `adata_metabolite.X`: float intensities, ideally already on a sane scale (post-PQN or post-log). The CCA and Spearman wrappers handle log1p; MMvec uses raw values directly.
- For `simulate_paired`: returns `(adata_mb, adata_mt, truth_df)` where `truth_df` columns are `microbe`, `metabolite`, indicating planted pairs.
- For Franzosa: writes CSVs under `data_dir`; first-run download is ~50 MB; subsequent runs read from cache.

## Minimal Execution Patterns

```python
import omicverse as ov
import numpy as np
import matplotlib.pyplot as plt
ov.plot_set()

# 1) Load real paired data (or use simulate_paired for testing)
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

# 3) Spearman triage
spear = ov.micro.paired_spearman(adata_mb_f, adata_mt_f, min_prevalence=0.1)
print(f'{(spear["fdr_bh"] < 0.05).sum()} pairs at FDR 0.05 / {len(spear)} tested')

# 4) CCA modes
cca = ov.micro.paired_cca(adata_mb_f, adata_mt_f, n_components=3)
print('canonical r:', [round(c, 3) for c in cca['canonical_correlations']])
ov.micro.plot_cca_scatter(cca, component=1)

# 5) MMvec
mmvec = ov.micro.MMvec(n_latent=4, epochs=600, val_frac=0.15,
                       patience=80, seed=0).fit(adata_mb_f, adata_mt_f)
print('best epoch:', mmvec.best_epoch_)
ov.micro.plot_mmvec_training(mmvec)
ov.micro.plot_cooccurrence(mmvec.cooccurrence(), top_n=15)
ov.micro.plot_embedding_biplot(mmvec, components=(0, 1), label_top=8)
print(mmvec.top_pairs(n=10))
```

```python
# Synthetic recovery — validate methods against planted pairs
ad_mb_sim, ad_mt_sim, truth = ov.micro.simulate_paired(n_pairs=5, seed=0)
res_sp = ov.micro.paired_spearman(ad_mb_sim, ad_mt_sim)
mmvec_sim = ov.micro.MMvec(n_latent=3, epochs=400, val_frac=0.1, seed=0).fit(
    ad_mb_sim, ad_mt_sim,
)
ov.micro.plot_paired_method_comparison(truth, spearman_df=res_sp, mmvec_model=mmvec_sim)
plt.tight_layout(); plt.show()
```

## Validation

- After loading: `adata_mb.obs_names.tolist() == adata_mt.obs_names.tolist()`. Mismatch causes `_check_paired` to raise.
- After Spearman: `len(spear) == adata_mb.n_vars * adata_mt.n_vars` (full pairwise). For 150 × 200 = 30k pairs, expect a few hundred at FDR 0.05.
- After CCA: `cca['canonical_correlations']` is a list of length `n_components`; first value should be larger than the rest. If first value < 0.5, the modalities don't share linear structure; CCA isn't the right method.
- After MMvec: `mmvec.best_epoch_ >= 0` (training found a checkpoint); `loss_history_` is monotonically decreasing on the train side; `val_loss_history_` shows a clear minimum at `best_epoch_`.
- If `val_loss_history_` plateaus while `loss_history_` keeps falling, you're overfitting — reduce `n_latent` or increase `l2`.
- For synthetic recovery: `plot_paired_method_comparison(truth, ...)` should show **planted pairs ranked higher than non-planted** under at least one method. If all three fail, your synthetic effect was buried below the noise floor — increase `effect_range` in `simulate_paired`.
- Cross-method agreement on real data: the top-10 by Spearman `|r|` and top-10 from `mmvec.top_pairs(10)` should overlap by at least 30 %. Zero overlap is a red flag (one of the methods misconfigured).

## Resource Map

- See [`reference.md`](reference.md) for compact copy-paste snippets.
- See [`references/source-grounding.md`](references/source-grounding.md) for verified `paired_*` / `MMvec.*` signatures and the docstring-fill log (`MMvec.fit / cooccurrence / conditional_probabilities / top_pairs / plot_mmvec_training / plot_embedding_biplot` — all backfilled by this skill's grounding pass).
- For microbe-side AnnData ingest, see `omicverse-microbiome-16s-amplicon-dada2`.
- For metabolite-side AnnData ingest + preprocessing, see `omicverse-bulk-metabol-preprocessing`.
- For running MOFA+ on metabol + RNA-seq (or three or more views), see `omicverse-bulk-metabol-pathway-multifactor`.

## Examples
- "Load Franzosa 2019 IBD paired data, filter to top 150 microbes / top 200 annotated metabolites, run Spearman + CCA + MMvec, and list the top-10 pairs from each method."
- "Validate MMvec on synthetic paired data with 5 planted producer-pairs (`simulate_paired(n_pairs=5)`); confirm planted pairs rank higher than random under MMvec and Spearman."
- "Plot the joint microbe + metabolite embedding biplot from a fitted MMvec at latent components (0, 1) with the top 8 entities labelled."
- "For a chosen microbe, return the conditional probability over metabolites under MMvec — i.e. `mmvec.conditional_probabilities().loc['<microbe_name>']`."

## References
- Tutorial notebook: [`t_micro_metabol_paired.ipynb`](https://omicverse.readthedocs.io/en/latest/Tutorials-Multi-Omics/micro-meta/t_micro_metabol_paired/) — Franzosa 2019 IBD walkthrough + synthetic validation.
- Live API verified — see [`references/source-grounding.md`](references/source-grounding.md).
