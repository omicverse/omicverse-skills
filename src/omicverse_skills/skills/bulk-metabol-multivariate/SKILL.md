---
name: omicverse-bulk-metabol-multivariate
description: Multivariate discrimination and biomarker selection on a preprocessed metabolomics AnnData. Use when running PLS-DA, OPLS-DA, VIP / S-plot inspection, per-metabolite ROC AUC with bootstrap confidence intervals, or building a multi-metabolite biomarker panel with nested CV and a permutation null. Assumes the input has already been imputed, PQN-normalized, and log-then-Pareto-transformed.
---

# OmicVerse Bulk Metabolomics — Multivariate & Biomarker

## Goal

Given a preprocessed metabolomics `AnnData` (output of the `omicverse-bulk-metabol-preprocessing` skill), discriminate two groups with PLS-DA / OPLS-DA, rank features by VIP, inspect the OPLS-DA S-plot, then build and validate a biomarker panel via per-feature ROC AUC and nested-CV multi-metabolite classification with a permutation null. The output is a ranked feature list, model-quality metrics (R²X / R²Y / Q² / mean AUC / permutation p-value), and an actionable biomarker panel.

This skill is **discrimination + biomarker selection only**. Imputation / normalization / transformation belong in the preprocessing skill; pathway interpretation of the resulting hits belongs in the pathway-multifactor skill.

## Quick Workflow

1. Confirm the input `AnnData` has been **PQN-normalized + log-transformed + Pareto-scaled**. Multivariate models on un-Pareto-scaled data are dominated by a few high-intensity features.
2. Fit `ov.metabol.plsda(...)` first as a baseline; inspect `r2x`, `r2y`, `q2`. A positive `Q²` is the minimum bar for a model that beats mean prediction.
3. Fit `ov.metabol.opls_da(n_ortho=1)` to separate group-discriminating signal (predictive component) from within-group variation (orthogonal components). Compare `r2y` / `q2` against the PLS-DA baseline.
4. Pull the VIP table with `result.to_vip_table(adata.var_names)` and visualize the top features with `ov.metabol.vip_bar(...)`.
5. Run `ov.metabol.s_plot(...)` to see covariance vs. correlation against the predictive component — features in the corners are robust biomarker candidates.
6. For biomarker selection: get per-feature univariate AUC + bootstrap CI with `ov.metabol.roc_feature(...)`.
7. Build a multi-metabolite panel with `ov.metabol.biomarker_panel(...)` using nested CV; include a permutation null to test whether the panel beats random feature selection.
8. Plot the held-out ROC curve from `panel.outer_predictions` / `panel.outer_labels` for the report.

## Interface Summary

Multivariate fitting (both consume Pareto-scaled `adata.X`):
- `ov.metabol.plsda(adata, *, group_col='group', group_a=None, group_b=None, n_components=2, scale=False) → PLSDAResult` — Partial Least Squares Discriminant Analysis.
- `ov.metabol.opls_da(adata, *, group_col='group', group_a=None, group_b=None, n_ortho=1, scale=False, max_iter=500, tol=1e-08) → PLSDAResult` — Orthogonal PLS-DA. Sets `scores` (predictive) and `x_ortho_scores` (orthogonal) separately.

`PLSDAResult` attributes:
- `.scores` — `(n_samples, n_components)` predictive scores.
- `.x_ortho_scores` — `(n_samples, n_ortho)` orthogonal scores (OPLS-DA only).
- `.loadings` — `(n_features, n_components)`.
- `.vip` — `(n_features,)` Variable Importance in Projection.
- `.coef` — `(n_features,)` regression coefficient against the binary Y (used for VIP-bar colour).
- `.r2x`, `.r2y`, `.q2` — standard PLS-DA quality metrics.
- `.group_labels` — tuple `(group_a, group_b)` actually fitted.
- `.to_vip_table(var_names) → pd.DataFrame` — VIP-sorted DataFrame with `vip`, `coef`, ready for biomarker triage.

Plotting:
- `ov.metabol.vip_bar(result, var_names, *, top_n=15, ax=None, figsize=(5, 5))` — horizontal VIP bars; bar colour encodes the sign of the regression coefficient (red = up in `group_a`, blue = up in `group_b`); dashed line at VIP=1 (Wold's threshold).
- `ov.metabol.s_plot(result, adata, *, label_top_n=15, ax=None, figsize=(5.5, 4.5))` — covariance vs. correlation against the predictive component; corner features are stable, high-magnitude markers.

Per-feature ROC:
- `ov.metabol.roc_feature(adata, *, group_col, pos_group=None, neg_group=None, layer=None, ci=False, n_bootstrap=1000, seed=0) → pd.DataFrame` with `auc`, optional `auc_ci_lower` / `auc_ci_upper`, plus per-feature `pvalue` against `auc=0.5`.

Multi-metabolite panel (nested CV):
- `ov.metabol.biomarker_panel(adata, *, group_col, features=10, classifier='rf'|'lr'|'svm', cv_outer=5, cv_inner=3, pos_group=None, neg_group=None, n_permutations=0, layer=None, seed=0) → BiomarkerPanelResult`.

`BiomarkerPanelResult` attributes:
- `.features` — final-panel metabolite names.
- `.classifier` — `'rf'`, `'lr'`, or `'svm'`.
- `.cv_outer`, `.cv_inner` — fold counts.
- `.outer_aucs` — `(cv_outer,)` held-out AUC per outer fold.
- `.outer_predictions`, `.outer_labels` — out-of-fold prediction arrays for ROC plotting.
- `.mean_auc`, `.std_auc` — summary across outer folds.
- `.permutation_pvalue` — `nan` when `n_permutations=0`; otherwise the empirical p-value vs. the label-permuted null.
- `.feature_importance` — DataFrame ranking selected features.
- `.pos_group`, `.neg_group` — fitted labels.
- `.to_frame() → pd.DataFrame` — per-fold AUC for report tables.

## Boundary

**Inside scope:**
- PLS-DA and OPLS-DA fitting on already-Pareto-scaled `AnnData`.
- VIP ranking, S-plot, score-plot inspection.
- Univariate per-feature AUC with bootstrap CI.
- Nested-CV multi-metabolite panel with permutation null.

**Outside scope — separate skill:**
- Imputation, normalization, transformation (preprocessing skill).
- Three-or-more-group discrimination — use `anova` from the preprocessing skill or a different multi-class model (out of scope here).
- Pathway / MSEA enrichment of VIP hits (pathway-multifactor skill).
- Multi-factor designs (ASCA / mixed model / MEBA — multifactor skill).
- Running PLS-DA on raw / log-only data — refusing this is a feature, not a missing capability.

## Branch Selection

**PLS-DA vs. OPLS-DA**
- `plsda` for a quick discriminative baseline. Fastest to fit; a positive `Q²` from `plsda` already proves the groups separate.
- `opls_da(n_ortho=1)` when the cohort is heterogeneous (within-group variability competes with between-group). Cleaner VIP table and clearer S-plot because orthogonal variation is partitioned out. Tutorial `t_metabol_02` shows `r2y` jumping substantially under OPLS-DA on the cachexia data.
- Increase `n_ortho` only when you have evidence (e.g. one prominent confounder still present after first orthogonal component); each additional `n_ortho` consumes degrees of freedom and `Q²` will eventually drop.

**`scale=False` is the right default** when the data was already Pareto-scaled in preprocessing. Setting `scale=True` here double-scales and degrades both `R²Y` and `Q²`.

**`n_components` for PLS-DA**: 2 is enough for a binary discrimination view. Add a third only if the score plot suggests a structured residual.

**`features` argument of `biomarker_panel`**
- `features=int` (e.g. 10) — picks the top-`int` features inside the inner CV loop using univariate ranking; the same algorithm runs in every outer fold so feature selection is honest.
- `features=list[str]` — fixes the panel up-front; useful when you have a literature-defined biomarker list and want to validate its CV performance.

**`classifier` choice**
- `'lr'` — logistic regression. Default for small cohorts (n < 50). Coefficients are interpretable; pairs naturally with PQN+log+Pareto-scaled features.
- `'rf'` — random forest. Robust to feature scaling but interpretability comes from `feature_importance`, not coefficients.
- `'svm'` — RBF-kernel SVM. Non-linear; least interpretable; use only when the linear models don't separate.

**`n_permutations`**
- Set `>= 100` for any panel you'll report. With `n_permutations=0` you get a mean AUC but no calibrated false-positive control; a panel with `mean_auc=0.95` on `n=20` samples might still be no better than random.
- `n_permutations=1000` if compute is available — gives 0.001-resolution p-values.

**`ci=True` on `roc_feature`** when you want per-feature CIs in the report; it costs `n_bootstrap` bootstrap resamples per feature. For the cachexia 77-sample dataset, `n_bootstrap=500` runs in seconds.

## Input Contract

- `AnnData` must come from the preprocessing skill: `obs['group']` populated, `X` is **log + Pareto-scaled** for `plsda` / `opls_da`, **log-only** for `roc_feature` / `biomarker_panel`.
- `obs['group']` must take exactly two values for PLS-DA / OPLS-DA / `biomarker_panel`. For 3+ groups, slice or mask first.
- For `biomarker_panel` with nested CV at `cv_outer=5, cv_inner=3`, you need at minimum 15 samples per class — fewer than that and the inner folds become trivial.
- Don't pass sparse `adata.X`; densify first.

## Minimal Execution Patterns

```python
# Multivariate (assumes preprocessing already ran)
import omicverse as ov
ov.plot_set()

# X must be log + Pareto-scaled
adata = ov.metabol.read_metaboanalyst('human_cachexia.csv', group_col='Muscle loss')
adata = ov.metabol.normalize(adata, method='pqn')
adata = ov.metabol.transform(adata, method='log')
adata = ov.metabol.transform(adata, method='pareto', stash_raw=False)

pls = ov.metabol.plsda(adata, group_col='group', n_components=2, scale=False)
print(f'PLS-DA  R²X={pls.r2x:.3f}  R²Y={pls.r2y:.3f}  Q²={pls.q2:.3f}')

opls = ov.metabol.opls_da(adata, group_col='group', n_ortho=1, scale=False)
print(f'OPLS-DA R²X={opls.r2x:.3f}  R²Y={opls.r2y:.3f}  Q²={opls.q2:.3f}')

vip_df = opls.to_vip_table(adata.var_names)
ov.metabol.vip_bar(opls, adata.var_names, top_n=15)
ov.metabol.s_plot(opls, adata, label_top_n=10)
```

```python
# Biomarker selection (consumes log-only AnnData; do NOT Pareto-scale here)
adata_log = ov.metabol.read_metaboanalyst('human_cachexia.csv', group_col='Muscle loss')
adata_log = ov.metabol.impute(adata_log, method='qrilc', seed=0)
adata_log = ov.metabol.normalize(adata_log, method='pqn')
adata_log = ov.metabol.transform(adata_log, method='log')

# 1) per-feature AUC with bootstrap CI
auc = ov.metabol.roc_feature(
    adata_log, group_col='group',
    pos_group='cachexic', neg_group='control',
    ci=True, n_bootstrap=500, seed=0,
)

# 2) nested-CV multi-metabolite panel
panel = ov.metabol.biomarker_panel(
    adata_log, group_col='group',
    pos_group='cachexic', neg_group='control',
    features=10,
    classifier='lr',
    cv_outer=5, cv_inner=3,
    n_permutations=100,
    seed=0,
)
print(f'panel mean AUC: {panel.mean_auc:.3f} ± {panel.std_auc:.3f}')
print(f'permutation p : {panel.permutation_pvalue:.3f}')
print(panel.feature_importance.head())

# 3) ROC curve from out-of-fold predictions
from sklearn.metrics import roc_curve, auc as _auc
import matplotlib.pyplot as plt
fpr, tpr, _ = roc_curve(panel.outer_labels, panel.outer_predictions)
fig, ax = plt.subplots(figsize=(5, 5))
ax.plot(fpr, tpr, lw=2, label=f'nested CV AUC = {_auc(fpr, tpr):.3f}')
ax.plot([0, 1], [0, 1], '--', color='gray')
ax.set(xlabel='FPR', ylabel='TPR', aspect='equal'); ax.legend(loc='lower right')
plt.tight_layout(); plt.show()
```

## Validation

- After fitting: `pls.q2 > 0` is the absolute minimum bar — negative `Q²` means the model is worse than predicting the mean. Any reportable result needs `Q² > 0.1` and ideally a permutation-derived empirical p-value.
- After OPLS-DA: `opls.r2y > pls.r2y` should hold on heterogeneous cohorts; if they're equal, the orthogonal direction didn't capture confounder structure and OPLS-DA is overkill.
- After `biomarker_panel`: report the *outer*-fold AUC (`mean_auc`), never the inner-fold AUC (which is contaminated by feature selection).
- `permutation_pvalue` should be < 0.05 for a publication-grade panel; if it's higher, the apparent discrimination is plausibly random.
- Confirm `panel.outer_labels` lines up with `adata.obs['group']` (the panel uses out-of-fold predictions, so labels are not in the original sample order — use the result's own arrays).
- For VIP-bar / S-plot interpretation: a feature being VIP > 1 alone is not enough — pair it with non-zero `coef` sign and a non-zero correlation in S-plot before calling it a biomarker.
- If the cohort is < 15 per class, the nested CV is unstable; fall back to univariate `roc_feature` only and don't claim panel performance.

## Resource Map

- See [`reference.md`](reference.md) for copy-paste snippets.
- See [`references/source-grounding.md`](references/source-grounding.md) for the verified `_plsda` / `_biomarker` / `plotting` signatures and notes on interpreting `PLSDAResult` and `BiomarkerPanelResult` attributes against tutorial usage.
- For the preprocessing chain that produces the input `AnnData`, see the `omicverse-bulk-metabol-preprocessing` skill.
- For pathway-level interpretation of the VIP / panel features, see the `omicverse-bulk-metabol-pathway-multifactor` skill.

## Examples
- "Run PLS-DA and OPLS-DA on the Pareto-scaled cachexia data and report `R²Y`, `Q²` for both."
- "Show the top-15 VIP metabolites with sign-coloured bars and an OPLS-DA S-plot."
- "Build a 10-metabolite logistic-regression biomarker panel with 5×3 nested CV and a 100-iteration permutation null; plot the held-out ROC curve."
- "Get per-feature AUC with 500-iteration bootstrap confidence intervals for cachexic vs. control."

## References
- Tutorial notebooks:
  - [`t_metabol_02_multivariate.ipynb`](https://omicverse.readthedocs.io/en/latest/Tutorials-metabol/t_metabol_02_multivariate/) — PLS-DA, OPLS-DA, VIP, S-plot.
  - [`t_metabol_08_biomarker.ipynb`](https://omicverse.readthedocs.io/en/latest/Tutorials-metabol/t_metabol_08_biomarker/) — `roc_feature`, nested-CV panel, permutation null.
- Live API verified against `omicverse.metabol` — see [`references/source-grounding.md`](references/source-grounding.md).
