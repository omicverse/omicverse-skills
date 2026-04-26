---
name: omicverse-microbiome-da-comparison
description: Run all three differential-abundance methods (Wilcoxon, pyDESeq2, ANCOM-BC) on the same microbiome AnnData, compare their hit sets via 3-way Venn / overlap counts, and decide which to trust on a given cohort. Use when the user wants to benchmark DA methods on a 16S study, when picking between methods on a small or zero-inflated cohort, or when reporting consensus features that survive multiple tests.
---

# OmicVerse Microbiome — Differential-Abundance Method Comparison

## Goal

Take a preprocessed microbiome `AnnData` (output of the 16S amplicon skill — samples × ASVs with 7-rank SINTAX taxonomy in `var`) and run **all three DA methods exposed by `ov.micro.DA`** on the same two-group contrast — Wilcoxon (rank), pyDESeq2 (NB-GLM), ANCOM-BC (compositional). Compare their hit sets at a common FDR cutoff, report consensus / method-specific genera, and surface the biology around the three methods' different statistical assumptions so the user can pick (or report) the right one.

## Quick Workflow

1. Load the AnnData produced by `omicverse-microbiome-16s-amplicon-dada2`. Drop control / non-relevant groups; keep exactly two for the contrast.
2. Collapse to a chosen taxonomic rank (typically genus): `ov.micro.collapse_taxa(adata, rank='genus')`. DA at species/ASV level is noisy on small cohorts; genus is the canonical reporting rank for 16S.
3. **Wilcoxon**: `ov.micro.DA(adata_genus).wilcoxon(group_key, group_a, group_b, min_prevalence=0.1)`. Non-parametric; fastest; tests ranks of relative abundance.
4. **pyDESeq2**: `ov.micro.DA(adata_genus).deseq2(group_key, group_a, group_b, min_prevalence=0.1)`. NB-GLM on **raw counts**; uses RNA-seq-style shrinkage of small-count fold-changes.
5. **ANCOM-BC**: `ov.micro.DA(adata_genus).ancombc(group_key, min_prevalence=0.1, pseudocount=1.0)`. Compositional bias-corrected ANCOM; closest to the compositional ground truth.
6. Build sets of significant features at a common FDR cutoff (typically `0.05`). Watch for column-name differences across methods (`fdr_bh` for Wilcoxon / DESeq2; `q_value` or `fdr_bh` for ANCOM-BC).
7. Tabulate the 3-way Venn (Wilcoxon-only / DESeq2-only / ANCOM-BC-only / pairwise overlaps / all-three) and render with `matplotlib_venn.venn3` if installed.
8. Report two numbers in the writeup: (a) the *consensus* hit set (intersection across all three) for the strongest claim, and (b) the *Wilcoxon ∪ ANCOM-BC* set as a slightly more permissive convention if compositional correctness matters.

## Interface Summary

Same `ov.micro.DA` class, three different methods, all returning a `pd.DataFrame` indexed by the feature column with method-specific columns:

```python
ov.micro.DA(adata).wilcoxon(
    group_key, group_a=None, group_b=None,
    rank=None, relative=True, min_prevalence=0.1,
) -> pd.DataFrame
# columns: feature, log2fc, pvalue, fdr_bh, mean_a, mean_b, prevalence_a, prevalence_b
# convention: log2fc > 0 ⇒ higher in group_a

ov.micro.DA(adata).deseq2(
    group_key, group_a=None, group_b=None,
    rank=None, min_prevalence=0.1, alpha=0.05,
) -> pd.DataFrame
# columns: feature, baseMean, log2fc, lfcSE, stat, pvalue, fdr_bh
# convention: log2fc > 0 ⇒ higher in group_a; uses pyDESeq2's MLE-shrunk LFC.

ov.micro.DA(adata).ancombc(
    group_key, rank=None,
    min_prevalence=0.1, pseudocount=1.0,
) -> pd.DataFrame
# columns: feature, log2fc, std_error, w_stat, pvalue, q_value (or fdr_bh)
# convention: log2fc > 0 ⇒ higher in group_a (after compositional bias correction).
```

## Boundary

**Inside scope:**
- Running all three DA methods on the same AnnData and the same two-group contrast.
- Set-arithmetic comparison of hit sets at a common FDR.
- Method-choice guidance based on cohort size / zero inflation / compositional concerns.
- Documenting consensus and method-specific features.

**Outside scope:**
- Building / preprocessing the input AnnData — see `omicverse-microbiome-16s-amplicon-dada2`.
- Cross-cohort meta-analysis (combining DA results from multiple studies) — see `omicverse-microbiome-meta-analysis`.
- Three-or-more-group DA — `ov.micro.DA` is two-group; use Kruskal-Wallis externally.
- Repeated-measures / paired DA — out of scope for `ov.micro` (use mixed models in `ov.metabol` or external `nlme` etc.).
- Phylogenetically-aware DA (e.g. PERMANOVA on UniFrac distances) — see phylogeny skill.

## Branch Selection

**Wilcoxon (default for moderate cohorts)**
- Pros: non-parametric (no NB / log-normal assumption); fast (sub-second on a 22-sample × 250-genus cohort); robust to outliers.
- Cons: weak power on small cohorts (n<10/group); doesn't model compositional bias; relative-abundance scaling makes interpretation of log2fc less rigorous than DESeq2.
- Use when: n>=10/group, low-to-moderate sparsity, no strong compositional concern, you need fast / portable.

**pyDESeq2 (RNA-seq-style NB-GLM)**
- Pros: well-calibrated p-values on raw counts; LFC shrinkage stabilises small-count fold-changes; mature implementation with extensive vignettes.
- Cons: NB assumption fails on highly zero-inflated data (microbiome can be 80%+ zeros); doesn't correct for compositional bias; slower than Wilcoxon.
- Use when: n is moderate, sparsity is moderate, you want NB-style shrinkage, reviewers expect DESeq2-style methodology.

**ANCOM-BC (compositional bias-corrected)**
- Pros: explicit compositional correction (the only method here that doesn't ignore the simplex); recovers true effect direction on highly biased data; bias-corrected per-feature log-ratio model is closest to the ground-truth.
- Cons: slowest of the three; requires `skbio>=0.7.1`; can be conservative on small cohorts; the bias-correction adds a per-feature constant which shifts log2fc relative to DESeq2.
- Use when: compositional correctness matters (always for true microbiome relative-abundance reporting), n is moderate-to-large, downstream analysis uses log-ratio interpretation.

**Consensus reporting strategy**
- For the *strongest* claim: intersection across all three (rare but unambiguous; survives method assumptions).
- For a *defensible* claim: Wilcoxon ∩ ANCOM-BC (combines a model-free test with the compositional-aware test; bypasses DESeq2's NB assumption).
- For *exploratory* hits: union; flag method-specific hits explicitly so the reader knows the assumption that drove them.
- Always report cohort size + sparsity + chosen FDR — DA results without those numbers are uninterpretable.

**Column-name pitfalls**
- Wilcoxon / pyDESeq2 use `fdr_bh`; ANCOM-BC may use `q_value` or `fdr_bh` depending on `skbio` version. Pattern in the tutorial: `sig_col = 'q_value' if 'q_value' in ab.columns else 'fdr_bh'`.
- All three methods may also expose a raw `pvalue` column. Don't confuse pvalue and FDR when filtering — pvalue<0.05 is *not* FDR<0.05.

## Input Contract

- An `AnnData` from the 16S amplicon skill; `obs[group_key]` is a categorical with at least two values; the two-group slice should have `n>=5/group` for any of these tests to behave reasonably.
- `adata.X` is **integer counts** (DESeq2 requires raw counts; Wilcoxon and ANCOM-BC handle either, but raw counts are the canonical input).
- For ANCOM-BC: `pip install skbio>=0.7.1` (function raises `ImportError` if missing).
- For pyDESeq2: `pip install pydeseq2`.

## Minimal Execution Patterns

```python
import omicverse as ov
import anndata as ad
import matplotlib.pyplot as plt

ov.plot_set()

# 1) Load AnnData from the 16S amplicon skill, restrict to two groups.
adata = ad.read_h5ad('mothur_sop_16s.h5ad')
adata = adata[adata.obs['group'].isin(['Early', 'Late'])].copy()
print(adata.obs['group'].value_counts().to_dict())

# 2) Collapse to genus.
adata_genus = ov.micro.collapse_taxa(adata, rank='genus')
print('samples × genera:', adata_genus.shape)

# 3) Run all three DA methods on the same contrast.
wx = ov.micro.DA(adata_genus).wilcoxon(
    group_key='group', group_a='Early', group_b='Late', min_prevalence=0.1,
)
ds = ov.micro.DA(adata_genus).deseq2(
    group_key='group', group_a='Early', group_b='Late', min_prevalence=0.1,
)
ab = ov.micro.DA(adata_genus).ancombc(
    group_key='group', min_prevalence=0.1,
)

print(f'  Wilcoxon : {(wx["fdr_bh"] < 0.05).sum():3d} / {len(wx):3d} genera at FDR 0.05')
print(f'  DESeq2   : {(ds["fdr_bh"] < 0.05).sum():3d} / {len(ds):3d}')
sig_col_ab = 'q_value' if 'q_value' in ab.columns else 'fdr_bh'
print(f'  ANCOM-BC : {(ab[sig_col_ab] < 0.05).sum():3d} / {len(ab):3d}')

# 4) Build sets and tabulate the 3-way Venn.
sig_wx = set(wx.loc[wx['fdr_bh']  < 0.05, 'feature'])
sig_ds = set(ds.loc[ds['fdr_bh']  < 0.05, 'feature'])
sig_ab = set(ab.loc[ab[sig_col_ab] < 0.05, 'feature'])

print('Wilcoxon only :', len(sig_wx - sig_ds - sig_ab))
print('DESeq2 only   :', len(sig_ds - sig_wx - sig_ab))
print('ANCOM-BC only :', len(sig_ab - sig_wx - sig_ds))
print('All three     :', len(sig_wx & sig_ds & sig_ab))

# 5) Render Venn (optional dependency)
try:
    from matplotlib_venn import venn3
    fig, ax = plt.subplots(figsize=(5, 5))
    venn3([sig_wx, sig_ds, sig_ab],
          set_labels=('Wilcoxon', 'DESeq2', 'ANCOM-BC'), ax=ax)
    ax.set_title('Genera significant at FDR/q < 0.05')
    plt.show()
except ImportError:
    pass
```

## Validation

- For each method: report `n_tested` (pre-min_prevalence filter), `n_significant` at FDR/q < 0.05, and the percentage. With a 22-sample mothur SOP demo, expect ~5–25 genus-level hits per method.
- Hit-set agreement: at least the top-3 by absolute log2fc per method usually agree; if Wilcoxon and ANCOM-BC have *zero* overlap, you have a problem (cohort too small, or the contrast doesn't actually exist in the data).
- pyDESeq2 hits that don't appear in Wilcoxon are often spurious driven by NB shrinkage on rare features; check the raw counts for those genera and consider tighter `min_prevalence`.
- ANCOM-BC's `log2fc` is bias-corrected — it may differ in *magnitude* from DESeq2's by a small constant offset (per the compositional correction), but the **sign** must agree on consensus hits. Sign disagreement on a consensus genus is a red flag.
- Reporting: always disclose cohort size, the FDR cutoff, and which `min_prevalence` was used. A "significant" feature at `min_prevalence=0.0` on n=5/group is uninterpretable.

## Resource Map

- See [`reference.md`](reference.md) for compact copy-paste snippets.
- See [`references/source-grounding.md`](references/source-grounding.md) for verified `DA.wilcoxon` / `DA.deseq2` / `DA.ancombc` signatures and column-name conventions across versions.
- For the AnnData ingest that produces the input, see `omicverse-microbiome-16s-amplicon-dada2`.
- For meta-analysis combining DA results across multiple cohorts, see `omicverse-microbiome-meta-analysis`.

## Examples
- "Run Wilcoxon, pyDESeq2, and ANCOM-BC on Early-vs-Late at the genus level and report the 3-way Venn at FDR 0.05."
- "Pick a single DA method for an n=12/group cohort with high zero-inflation — give the rationale."
- "Build the consensus hit set (intersection) and the Wilcoxon ∩ ANCOM-BC set for a publication report."
- "Diagnose why a DESeq2 hit doesn't appear in Wilcoxon — likely an NB-shrinkage artefact on a low-prevalence feature."

## References
- Tutorial notebook: [`t_16s_da_comparison.ipynb`](https://omicverse.readthedocs.io/en/latest/Tutorials-microbiome/t_16s_da_comparison/) — three-way DA Venn on the mothur SOP cohort.
- Live API verified — see [`references/source-grounding.md`](references/source-grounding.md).
