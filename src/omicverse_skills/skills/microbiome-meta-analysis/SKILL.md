---
name: omicverse-microbiome-meta-analysis
description: Combine multiple per-study microbiome AnnDatas into a single cross-cohort table and run inverse-variance / random-effects meta-analysis on differential abundance. Use when you have 16S studies from multiple cohorts that need joint analysis, when you need a combined log2 fold-change with Cochran's I² heterogeneity, or when you want to find features whose effect replicates across cohorts (vs. cohort-specific signals).
---

# OmicVerse Microbiome — Cross-cohort Meta-analysis

## Goal

Take per-study microbiome `AnnData` objects (one per cohort, each from the 16S amplicon skill), stitch them into a single cross-cohort AnnData with `combine_studies`, then run **per-study DA → cross-study meta-analysis** with `meta_da`. Output is a single feature-level table with combined effect size, combined SE, p-value, BH-FDR, **per-study count**, and Cochran's I² heterogeneity — distinguishing features whose effect *replicates* across cohorts from cohort-specific signals.

The skill assumes the meta question is **"does this contrast hold across cohorts?"**, not "is there an effect in any cohort". The latter is just running per-study DA and unioning hits — call that out as the wrong tool.

## Quick Workflow

1. Load (or build) per-cohort `AnnData` — same conventions as the 16S amplicon skill (samples × ASVs, 7-rank SINTAX in `var`, `obs[group_key]` for the contrast).
2. **Combine**: `ov.micro.combine_studies([a, b, c], study_names=[...], rank='genus', study_key='study') → AnnData` — collapses ASVs to a common rank (genus by default), unions features across studies (zero-fills cohort-specific features), tags `obs['study']`.
3. **Diagnostic Beta + PCoA** on the combined object: confirms whether study effect dominates the disease effect — a useful sanity check before running meta_da. If PC1 colour-maps cleanly by `study`, the cohort effect is large and meta_da's random-effects pooling is doing real work.
4. **Per-study top hits** (optional, for sanity): run `DA(a_g).deseq2(...)` on each study independently, take top-5 by p-value; check whether the planted / known-true features appear in every cohort.
5. **Meta-analysis**: `ov.micro.meta_da([a, b, c], study_names=[...], group_key, group_a, group_b, method='deseq2'|'wilcoxon'|'ancombc', rank='genus', min_prevalence=0.1, combine='random_effects'|'fixed_effects') → pd.DataFrame` with `feature`, `combined_lfc`, `combined_se`, `z`, `p_value`, `fdr_bh`, `n_studies`, `I2`.
6. **Plot the result**: forest-style barh of top features (combined LFC ± 1.96·SE), colour-coding planted vs. cohort-specific vs. shared-noise features.
7. **Heterogeneity scatter**: combined LFC × I² — high-I² features (>0.75) are reproducibly different in *direction* between cohorts (suspicious); low-I² features (<0.25) are coherent. Use this to flag which "significant" hits actually replicate.
8. (Optional) **Method × meta benchmark**: re-run meta_da with `method` cycled over `'wilcoxon' / 'deseq2' / 'ancombc'`, compare hit counts and the rank of known-true features across methods.

## Interface Summary

```python
ov.micro.combine_studies(
    studies: Sequence[AnnData],
    study_names: Optional[Sequence[str]] = None,
    rank: Optional[str] = 'genus',         # collapse to this rank first; None keeps ASVs
    study_key: str = 'study',              # column name in combined obs
    min_prevalence: float = 0.0,           # additional filter applied after concatenation
) -> AnnData

ov.micro.meta_da(
    studies: Sequence[AnnData],
    group_key: str,
    group_a: Optional[str] = None,
    group_b: Optional[str] = None,
    method: str = 'deseq2',                # 'wilcoxon' | 'deseq2' | 'ancombc'
    rank: Optional[str] = 'genus',
    min_prevalence: float = 0.1,
    combine: str = 'random_effects',       # or 'fixed_effects'
    study_names: Optional[Sequence[str]] = None,
    **method_kwargs,
) -> pd.DataFrame
# columns: feature, combined_lfc, combined_se, z, p_value, fdr_bh,
#          n_studies, I2, plus per-study lfc / se / pvalue columns
```

## Boundary

**Inside scope:**
- Combining 2+ per-study AnnDatas at any taxonomic rank.
- Per-study DA + cross-study inverse-variance pooling with random- or fixed-effects model.
- Cochran's I² heterogeneity diagnostics.
- Diagnostic combined-PCoA to spot study effects.
- Method-comparison sweep for meta-analysis (wilcoxon / deseq2 / ancombc).

**Outside scope — separate skill:**
- Building the per-study AnnData from FASTQs — see `omicverse-microbiome-16s-amplicon-dada2`.
- Single-cohort DA method comparison — see `omicverse-microbiome-da-comparison`.
- Phylogenetic harmonisation across studies (different reference trees) — see `omicverse-microbiome-phylogeny` for the per-study tree, but no cross-study tree alignment is implemented.
- Batch correction in the canonical sense (ComBat / Harmony) — `combine_studies` only zero-fills and tags `obs['study']`. For batch correction across cohorts, consider modelling `study` as a covariate in a per-feature mixed model (out of scope for `ov.micro`).
- Strain-level meta-analysis — current implementation collapses to a taxonomic rank; species-level cross-cohort matching needs uniform DBs (out of scope).

## Branch Selection

**`combine_studies(rank=...)`**
- `rank='genus'` (default) — robust across cohorts, the canonical reporting rank for 16S meta-analyses (Nearing 2022, Pasolli 2017). Use this unless you have a specific reason.
- `rank='family'` — when even genera don't align well across DBs.
- `rank=None` — ASV-level; only useful if every cohort was processed against the same DB version with the same denoiser; rare in practice.

**`meta_da(method=...)`**
- `'deseq2'` — default. Per-study NB-GLM, then inverse-variance pooling on the per-study LFCs. Best calibrated when `n>=10/group` per cohort.
- `'wilcoxon'` — non-parametric per-study, then pool. Use when one cohort has very small `n` where NB assumptions break.
- `'ancombc'` — compositional bias-corrected per-study, then pool. The most rigorous on compositional data but slowest; flag for reportable results.

**`combine='random_effects'` vs `'fixed_effects'`**
- `random_effects` (DerSimonian-Laird, default): models a between-study variance component τ² in addition to per-study sampling variance. Conservative when studies disagree (high I²); reduces to fixed-effects when studies agree.
- `fixed_effects`: assumes one true effect across studies; tighter CIs but mis-leads when studies actually differ. Use only when between-study heterogeneity is provably zero (rare for 16S).
- Default to random-effects for 16S meta-analyses — different cohorts have different recruitment, geography, sample handling.

**Cochran's I² interpretation**
- I² < 0.25: low heterogeneity — combined effect is well-supported.
- 0.25 ≤ I² < 0.5: moderate; report combined effect but flag.
- 0.5 ≤ I² < 0.75: substantial; combined effect is suspect — does the contrast actually hold across cohorts?
- I² >= 0.75: high; almost certainly a method or covariate problem (or genuine biology that varies across cohorts) — interpret per-study, don't pool.

**`min_prevalence`** in meta_da is applied per-study (not after combination), so a feature must clear `min_prevalence` in **at least one** study to be tested. Setting it lower (0.05) tests more rare features but inflates FDR-correction burden.

## Input Contract

- Per-study AnnDatas all have:
  - `obs[group_key]` populated with categorical values; the two contrast levels (`group_a`, `group_b`) appear in **every** study (otherwise that study contributes nothing to the contrast).
  - `var[rank]` for whatever taxonomic rank `combine_studies` will collapse to. Empty strings (unassigned features) are dropped during collapse.
  - Integer counts in `.X` (DESeq2 mode) or counts that survive `relative=True` rescaling (Wilcoxon mode).
- Sample counts per study can be unequal — the random-effects model handles this; just don't have any cohort with `n<5/group`.
- Different per-study DBs are **OK** at the genus level (genus names are largely DB-portable). Below genus, alignment breaks.

## Minimal Execution Patterns

```python
import omicverse as ov
import anndata as ad
import matplotlib.pyplot as plt
import numpy as np

ov.plot_set()

# 1) Load per-study AnnDatas
adata_A = ad.read_h5ad('cohort_A_genus.h5ad')
adata_B = ad.read_h5ad('cohort_B_genus.h5ad')
adata_C = ad.read_h5ad('cohort_C_genus.h5ad')

# 2) Combine — zero-fills cohort-specific features, tags obs['study']
adata_all = ov.micro.combine_studies(
    [adata_A, adata_B, adata_C],
    study_names=['cohort_A', 'cohort_B', 'cohort_C'],
    rank='genus',
)
print('combined:', adata_all.shape)
print(adata_all.obs['study'].value_counts())

# 3) Diagnostic PCoA on the combined object
ov.micro.Beta(adata_all).run(metric='braycurtis', rarefy=False)
coords = np.asarray(ov.micro.Ordinate(adata_all, dist_key='braycurtis').pcoa(n=2))

# 4) Meta-analysis
meta = ov.micro.meta_da(
    [adata_A, adata_B, adata_C],
    study_names=['cohort_A', 'cohort_B', 'cohort_C'],
    group_key='disease', group_a='CTRL', group_b='CASE',
    method='deseq2',
    rank='genus',
    min_prevalence=0.1,
    combine='random_effects',
)

print(f'features tested cross-cohort: {len(meta)}')
print(f'significant at FDR 0.05    : {(meta["fdr_bh"] < 0.05).sum()}')
print(meta[['feature', 'combined_lfc', 'combined_se', 'z',
            'p_value', 'fdr_bh', 'n_studies', 'I2']].head(10))
```

```python
# Forest-style top hits with 95 % CI (1.96 × combined_se)
top = meta.head(12).copy().iloc[::-1]
ci = 1.96 * top['combined_se'].values

fig, ax = plt.subplots(figsize=(7, 5))
ax.barh(top['feature'], top['combined_lfc'], xerr=ci,
        edgecolor='black', linewidth=0.4)
ax.axvline(0, color='k', lw=0.7)
ax.set_xlabel('combined log2 fold-change (CASE / CTRL)')
plt.tight_layout(); plt.show()
```

```python
# Heterogeneity scatter — flag high-I² hits
sig = meta[meta['fdr_bh'] < 0.05].copy()
fig, ax = plt.subplots(figsize=(7, 4))
ax.scatter(sig['combined_lfc'], sig['I2'], s=50, edgecolors='k', linewidths=0.4)
ax.axhline(0.25, color='grey', lw=0.5, ls='--')
ax.axhline(0.75, color='grey', lw=0.5, ls='--')
ax.set_xlabel('combined log2FC'); ax.set_ylabel("Cochran's I²")
plt.tight_layout(); plt.show()
```

```python
# Per-method comparison — which method recovers planted signals first?
results = {}
for method in ['wilcoxon', 'deseq2', 'ancombc']:
    results[method] = ov.micro.meta_da(
        [adata_A, adata_B, adata_C],
        study_names=['cohort_A', 'cohort_B', 'cohort_C'],
        group_key='disease', group_a='CTRL', group_b='CASE',
        method=method, rank='genus', min_prevalence=0.1,
    )
import pandas as pd
summary = pd.DataFrame({
    'n_tested':       {m: len(r) for m, r in results.items()},
    'n_sig_fdr_0.05': {m: int((r['fdr_bh'] < 0.05).sum()) for m, r in results.items()},
})
print(summary)
```

## Validation

- After `combine_studies`: `adata_all.obs['study'].value_counts()` matches per-study sample counts; combined feature count is `<= sum(per_study_features)` (zero-fills shouldn't lose features).
- Diagnostic PCoA: PC1 colour-mapped by `study` should look at least somewhat banded — if not, you have miraculously homogeneous data, but check that `combine_studies` collapsed correctly (rank actually shared).
- `meta_da` returns `n_studies` per feature — features that only show up in 1–2 studies have weaker pooled inference. Flag rows where `n_studies < len(studies)`.
- I² above 0.75 on a "significant" feature: do **not** report as a meta-finding. Interpret per-study or note explicitly that the effect is heterogeneous.
- For known-truth validation (planted features in synthetic data, or housekeeping features in real data): they should appear at the top of the per-study DA AND survive meta-analysis with low I² and small p-value. If a planted feature has high I² in meta_da, the synthetic noise wasn't realistic enough.
- Cross-method comparison: for at least 2/3 methods, the planted features should be in the top 3 hits. Methods that completely miss planted features signal a configuration error (wrong contrast, wrong `min_prevalence`).

## Resource Map

- See [`reference.md`](reference.md) for compact copy-paste snippets per stage.
- See [`references/source-grounding.md`](references/source-grounding.md) for verified `combine_studies` / `meta_da` signatures and the inverse-variance + random-effects implementation details.
- For per-study AnnData ingest, see `omicverse-microbiome-16s-amplicon-dada2`.
- For single-cohort DA method comparison, see `omicverse-microbiome-da-comparison`.

## Examples
- "Combine three cohorts at the genus level, run a random-effects DESeq2 meta-analysis, and report the top-12 features with 95 % CIs and I² heterogeneity."
- "Compare Wilcoxon / DESeq2 / ANCOM-BC at the meta level — which recovers planted features first?"
- "Diagnose whether a 'significant' meta-hit at FDR 0.05 has acceptable cross-cohort agreement (I²)."
- "Run a fixed-effects meta-analysis only when all three studies agree (low I²)."

## References
- Tutorial notebook: [`t_16s_meta_analysis.ipynb`](https://omicverse.readthedocs.io/en/latest/Tutorials-microbiome/t_16s_meta_analysis/) — synthetic 3-cohort meta-analysis with planted vs. cohort-specific features.
- Live API verified — see [`references/source-grounding.md`](references/source-grounding.md).
