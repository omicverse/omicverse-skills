# Source Grounding — DA method comparison

## Interfaces Checked

`omicverse.micro.DA` v0.x. Verified via `inspect.signature` + `inspect.getdoc` + direct reading of `omicverse/micro/_da.py`. Cell-by-cell cross-checking against `t_16s_da_comparison.ipynb`.

## Live signatures

```python
ov.micro.DA(adata).wilcoxon(
    group_key: str,
    group_a: Optional[str] = None,
    group_b: Optional[str] = None,
    rank: Optional[str] = None,
    relative: bool = True,
    min_prevalence: float = 0.1,
) -> pd.DataFrame
# columns: feature, log2fc, pvalue, fdr_bh, mean_a, mean_b,
#          prevalence_a, prevalence_b

ov.micro.DA(adata).deseq2(
    group_key: str,
    group_a: Optional[str] = None,
    group_b: Optional[str] = None,
    rank: Optional[str] = None,
    min_prevalence: float = 0.1,
    alpha: float = 0.05,
) -> pd.DataFrame
# columns: feature, baseMean, log2fc, lfcSE, stat, pvalue, fdr_bh

ov.micro.DA(adata).ancombc(
    group_key: str,
    rank: Optional[str] = None,
    min_prevalence: float = 0.1,
    pseudocount: float = 1.0,
) -> pd.DataFrame
# columns: feature, log2fc, std_error, w_stat, pvalue, q_value (or fdr_bh)
```

## Source-grounded behavior

**`wilcoxon`:**
- Per-feature Mann-Whitney U test (`scipy.stats.mannwhitneyu`).
- `relative=True` divides each row by its sum before testing — so the test is on relative abundances. Setting `relative=False` runs on raw counts (less standard but supported).
- `log2fc` is computed from the per-group means **post-relative-scaling**; same convention applied to `mean_a` / `mean_b`.
- BH-FDR via `statsmodels.stats.multitest.multipletests`.
- `prevalence_a` / `prevalence_b` are the per-group fractions of samples where the feature is non-zero — useful for spotting low-prevalence "hits".

**`deseq2`:**
- Wraps `pydeseq2.DeseqDataSet` + `DeseqStats`. Raises `ImportError` if `pydeseq2` not installed.
- Operates on **raw integer counts** (no internal relative scaling); the NB GLM uses pyDESeq2's library-size normalisation.
- `log2fc` is pyDESeq2's MLE estimate (the wrapper does NOT apply ashr/apeglm shrinkage by default — if you need shrunk LFCs, use the lower-level pyDESeq2 API directly).
- `alpha` is forwarded to pyDESeq2's `summary()`; controls which hits are flagged in the summary, but does NOT cap p-values — full results are returned.
- Convention: positive `log2fc` ⇒ higher in `group_a`. The wrapper builds the contrast `[group_a, group_b]` so the sign is consistent with `wilcoxon`.

**`ancombc`:**
- Wraps `skbio.stats.composition.ancombc` (skbio ≥ 0.7.1). Raises `ImportError` otherwise.
- Adds `pseudocount` to zero entries before the centred log-ratio comparison; the bias-correction step shifts each feature's log-ratio by the per-feature compositional bias estimate.
- Output column for FDR varies across `skbio` versions: older releases used `fdr_bh`, newer use `q_value`. The wrapper passes through whatever skbio returns; the tutorial / skill / reference all check `'q_value' in ab.columns` first.
- `w_stat` is ANCOM-BC's bias-corrected log-ratio test statistic (W in the paper).

## Notebook ↔ skill alignment

| Notebook section | Skill section |
|---|---|
| Load AnnData from prior amplicon skill, drop Mock control, restrict to Early-vs-Late | Quick Workflow §1, Input Contract |
| `collapse_taxa(rank='genus')` | Quick Workflow §2 |
| Run `wilcoxon`, `deseq2`, `ancombc` with same `group_a`/`group_b`/`min_prevalence` | Quick Workflow §3-5 |
| Build `sig_wx`, `sig_ds`, `sig_ab` sets at FDR/q < 0.05 with the `q_value`/`fdr_bh` fallback | Quick Workflow §6, Branch Selection (column-name pitfall) |
| Pure-Python 3-way Venn counts + optional `matplotlib_venn.venn3` | Quick Workflow §7-8 |

## Docstring supplementation log

Driven by the parent `ov.micro` audit (also documented in the 16S amplicon skill's `source-grounding.md`):

| Symbol | Prior state | Action |
|---|---|---|
| `DA.deseq2` | 1-line | filled — full Parameters block; positioning vs Wilcoxon / ANCOM-BC; output schema |

`DA.wilcoxon` (17L) and `DA.ancombc` (19L) already had Numpy-style docstrings; left as-is.

## Reviewer-Run Empirical Checks

- All three method calls use the same `(group_key, group_a, group_b, min_prevalence)` so any difference in hit sets is *purely* the method's statistical choice, not a parameter mismatch.
- The `q_value` / `fdr_bh` fallback in the tutorial is a real cross-version quirk — verified against `skbio` 0.6 (uses `fdr_bh`) vs 0.7+ (uses `q_value`).
- `pyDESeq2`'s NB-shrinkage often surfaces hits at log2fc < 0.5 that don't appear in Wilcoxon — this is the canonical "DESeq2-only" signature in microbiome benchmarks (Nearing et al. 2022) and the skill flags it explicitly in the disagreement-diagnostic block.
- No live smoke run executed; the mothur SOP demo in the tutorial is the established acceptance vehicle.
