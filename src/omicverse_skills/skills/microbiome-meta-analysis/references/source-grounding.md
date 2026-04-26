# Source Grounding — Cross-cohort meta-analysis

## Interfaces Checked

`omicverse.micro.combine_studies` and `omicverse.micro.meta_da` (`_meta.py`). Verified via `inspect.signature` + `inspect.getdoc` + direct reading of `omicverse/micro/_meta.py`. Cell-by-cell cross-checking against `t_16s_meta_analysis.ipynb`.

## Live signatures

```python
ov.micro.combine_studies(
    studies: Sequence[ad.AnnData],
    study_names: Optional[Sequence[str]] = None,
    rank: Optional[str] = 'genus',
    study_key: str = 'study',
    min_prevalence: float = 0.0,
) -> ad.AnnData

ov.micro.meta_da(
    studies: Sequence[ad.AnnData],
    group_key: str,
    group_a: Optional[str] = None,
    group_b: Optional[str] = None,
    method: str = 'deseq2',
    rank: Optional[str] = 'genus',
    min_prevalence: float = 0.1,
    combine: str = 'random_effects',
    study_names: Optional[Sequence[str]] = None,
    **method_kwargs,
) -> pd.DataFrame
```

## Source-grounded behavior

**`combine_studies`:**
- Per-study `collapse_taxa(rank=...)` first (when `rank` is given), so all studies talk about the same taxonomic units before unioning.
- Unions feature names across studies; for each study, features missing from that study get **zero-filled** in the resulting matrix. The combined sparse matrix carries the cohort-specific zero pattern in its non-trivial structure.
- `obs['study']` (or whatever `study_key` says) is added as a categorical column; `study_names` defaults to `['study_0', 'study_1', ...]` if not provided.
- `min_prevalence` here is a global filter applied AFTER combining; useful to drop ultra-rare features that only appear in one study at trace levels. Defaults to 0.0 (keep everything zero-filled) so the user can decide.
- Combined `var` columns: only the rank column survives the collapse; other taxonomy ranks are filled with empty strings (since collapse loses sub-rank resolution).

**`meta_da`:**
- Per-study DA: dispatches to `DA(adata_i).<method>(group_key, group_a, group_b, rank, min_prevalence, **method_kwargs)`. Method-specific kwargs forward through.
- Inverse-variance pooling: for each feature appearing in `>=2` studies, the per-study log-fold-changes are pooled with weights `w_i = 1 / SE_i^2` (fixed-effects) or `w_i = 1 / (SE_i^2 + τ²)` (random-effects, DerSimonian-Laird τ²).
- `n_studies` is the per-feature count of studies where the feature passed `min_prevalence` and the per-study DA was actually computed.
- `I2 = max(0, (Q - df) / Q)` where `Q` is Cochran's Q and `df = n_studies - 1`. Bounded `[0, 1]`. Conservative (truncates negative I²).
- The output DataFrame also carries per-study `lfc_<study_name>` and `se_<study_name>` columns when `study_names` is supplied — useful for downstream forest plots.
- BH-FDR via `statsmodels.stats.multitest.multipletests` over all features that returned a finite combined p-value.

**`combine='fixed_effects'`** assumes τ²=0 — tighter CIs but mis-leads when studies disagree. The default `random_effects` is the right choice for 16S meta-analyses where geography / handling / DB version vary across cohorts.

**`method` flow-through:**
- `method='deseq2'` requires `pip install pydeseq2` per study; the wrapper handles per-study failure gracefully (logs a warning and excludes that study from the per-feature pool).
- `method='ancombc'` requires `skbio>=0.7.1`. ANCOM-BC's bias-corrected log-ratio is what gets pooled (not raw log2fc), so combined values from ANCOM-BC may be slightly offset from DESeq2's even on the same data.
- `method='wilcoxon'` is fastest and most permissive on small per-cohort `n`.

## Notebook ↔ skill alignment

| Notebook section | Skill section |
|---|---|
| Synthetic 3-cohort `simulate_cohort` builder | Demoted — example construction; not in skill body |
| `combine_studies` with `rank='genus'` + per-study breakdown print | Quick Workflow §2; Validation |
| Combined Beta + PCoA showing study effect dominates | Quick Workflow §3 |
| Per-study `DA.deseq2` top-5 sanity check across cohorts | Quick Workflow §4 |
| `meta_da(method='deseq2', combine='random_effects')` + cross-cohort hits | Quick Workflow §5; Branch Selection |
| Forest-style barh of top features (LFC ± 1.96·SE) coloured by feature class | Quick Workflow §6; reference.md |
| Effect × I² heterogeneity scatter | Quick Workflow §7; reference.md |
| Per-method comparison sweep (wilcoxon / deseq2 / ancombc) | Quick Workflow §8; reference.md benchmark block |

## Docstring supplementation log

`combine_studies` (32L) and `meta_da` (47L) had complete Numpy-style docstrings prior to this skill — no supplementation needed. The synthetic cohort simulator in the tutorial is an example helper, not a public API.

## Reviewer-Run Empirical Checks

- `ov.micro.combine_studies` and `ov.micro.meta_da` both importable and runnable on the synthetic 3-cohort dataset constructed in the notebook (verified by reading the source paths; not by running).
- Random-effects τ² formula verified against DerSimonian-Laird's 1986 formulation in `_meta.py`.
- Output schema (`combined_lfc`, `combined_se`, `z`, `p_value`, `fdr_bh`, `n_studies`, `I2`) confirmed via the notebook's `meta[cols].head(10)` inspection.
- No live smoke run executed.
