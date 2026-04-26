---
name: omicverse-bulk-metabol-preprocessing
description: Run the canonical metabolomics preprocessing chain on an AnnData peak table — impute, normalize, transform — and apply LC-MS-specific drift / batch / sample QC corrections (drift_correct, SERRF, ComBat, sample_qc). Use when converting `t_metabol_01_intro` or `t_metabol_06_batch_correction` into a reusable skill, when a user has a MetaboAnalyst CSV / LC-MS peak table to clean before differential or multivariate analysis, or when choosing among PQN/TIC/median sample-normalization, log/Pareto feature-transformation, and qrilc/knn/half_min/zero imputation strategies.
---

# OmicVerse Bulk Metabolomics Preprocessing

## Goal

Take a metabolomics peak table (NMR or LC-MS; MetaboAnalyst CSV, generic wide CSV, or `m/z`/`RT`-coded LC-MS table) into an `AnnData` and run the canonical clean-up chain — load → impute → sample-normalize → feature-transform — plus optional MS-specific drift / batch / sample-outlier correction. The output is an analysis-ready `AnnData` with `obs['group']` populated, raw intensities preserved as a layer, and the matrix in the right scale for the downstream skill (univariate stats vs. multivariate model).

Stay focused on preprocessing. Univariate `differential` / `anova` is a thin add-on shown for completeness; multivariate (PLS-DA, OPLS-DA), pathway enrichment (MSEA), biomarker selection, DGCA, and multi-omics MOFA are separate downstream skills.

## Quick Workflow

1. Pick the loader by file format and pass `group_col` — every downstream test reads `adata.obs['group']`.
2. Impute missing values; method depends on the missingness mechanism (MAR → `knn`; MNAR / left-censored LC-MS → `qrilc`).
3. Normalize each sample row to remove dilution effects — PQN is canonical.
4. Transform each feature column for variance stabilization (`log`), and Pareto-scale on top **only** when a multivariate model will consume the matrix.
5. For LC-MS with QC pools and injection order: run `drift_correct` (LOESS-on-QC) → `serrf` (RF on QC) → `ov.bulk.batch_correction` (ComBat on residual batch shifts), in that order.
6. Run `sample_qc` on the corrected real samples to flag outliers; inspect manually before removing.
7. Validate that `adata.obs['group']` exists, no NaNs remain in `adata.X`, and `raw` matrix is preserved in a layer.

## Interface Summary

Loaders (each returns `AnnData` with `obs=samples`, `var=metabolites`, `obs['group']` populated):
- `ov.metabol.read_metaboanalyst(path, *, group_col, sample_col=None, transpose=False)` — MetaboAnalyst-format CSV.
- `ov.metabol.read_wide(path, *, sample_col, group_col, ...)` — generic wide CSV with explicit axes.
- `ov.metabol.read_lcms(path, *, feature_id_sep='/', sample_col=None, group_col=None, label_row=None, transpose=True)` — LC-MS peak table; parses `m/z` and `RT` from feature IDs into `var['mz']` / `var['rt']`.

Preprocessing chain (each is `AnnData`-in / `AnnData`-out and stashes the previous matrix into `adata.layers` when `stash_raw=True`):
- `ov.metabol.impute(adata, *, method='qrilc'|'knn'|'half_min'|'zero', missing_threshold=0.5, n_neighbors=5, q=0.01, seed=0)`.
- `ov.metabol.normalize(adata, *, method='pqn'|'tic'|'median'|'mstus', reference='median'|'mean', missing_threshold=0.5)`.
- `ov.metabol.transform(adata, *, method='log'|'glog'|'autoscale'|'pareto', pseudocount=1.0, stash_raw=True)`.

QC and MS-specific correction:
- `ov.metabol.cv_filter(adata, *, qc_mask=None, cv_threshold=0.30, across='qc'|'all')` — drop high-CV features.
- `ov.metabol.blank_filter(adata, *, blank_mask, ratio=3.0)` — drop features not above blank intensity.
- `ov.metabol.drift_correct(adata, *, injection_order, qc_mask, frac=0.5)` — LOESS regression on QC samples vs. injection order.
- `ov.metabol.serrf(adata, *, qc_col, qc_label='QC', batch_col=None, top_k=10, n_estimators=100, min_qc_samples=5, layer=None, seed=0)` — RF residual correction (Fan 2019). Adds `var['cv_qc_raw']` and `var['cv_qc_serrf']`.
- `ov.bulk.batch_correction(adata, *, batch_key, key_added='batch_correction')` — ComBat; corrected matrix lands at `adata.layers['batch_correction']`.
- `ov.metabol.sample_qc(adata, *, n_components=2, alpha=0.95, center=True, scale=True, layer=None) → pd.DataFrame` with `T2`, `DModX`, `T2_crit`, `DModX_crit`, `is_outlier`.
- `ov.metabol.sample_qc_plot(qc_df, *, ax=None, figsize=(5, 4), normal_color, outlier_color)` — Hotelling T² vs. DModX scatter with critical-value lines.

Univariate stats (kept here for completeness — they consume the log-transformed matrix, not the Pareto one):
- `ov.metabol.differential(adata, *, group_col='group', group_a, group_b, method='welch_t'|'mannwhitney'|'paired_t', layer=None, log_transformed=True) → pd.DataFrame` with `log2fc`, `pvalue`, `padj`, `mean_a`, `mean_b`.
- `ov.metabol.anova(adata, *, group_col='group', groups=None, method='welch_anova'|'anova'|'kruskal', layer=None) → pd.DataFrame`.
- `ov.metabol.volcano(deg, *, padj_thresh=0.05, log2fc_thresh=1.0, label_top_n=10, use_pvalue=False, clip_log2fc=None, ax=None, figsize=(5.5, 4.5))`.

Chained class API:
- `ov.metabol.pyMetabo(adata, random_state=0)` — keeps `self.raw` frozen, exposes `.cv_filter / .drift_correct / .blank_filter / .impute / .normalize / .transform / .differential / .plsda / .opls_da` (all chainable, all return `self`), plus `.deg_table`, `.plsda_result`, `.vip_table()`, `.significant_metabolites(padj_thresh, log2fc_thresh)`.

## Boundary

Keep this skill focused on the deterministic preprocessing chain plus MS-specific correction. **Inside scope:**
- Loading any of the three peak-table formats.
- Imputation, sample-normalization, feature-transformation in any order, with reproducible seeding.
- LC-MS drift / SERRF / ComBat / sample-QC stack.
- Demonstration that the `pyMetabo` chain produces identical numerics to the functional API (notebook asserts equality to machine precision).

**Outside scope — separate skill:**
- PLS-DA / OPLS-DA model fitting, VIP scoring, S-plot (multivariate skill).
- MSEA, mummichog, KEGG / LIPID MAPS / LION / HMDB ID mapping (pathway skill).
- Biomarker panel + nested CV (biomarker skill).
- DGCA differential correlation (correlation skill).
- ASCA / mixed model / MEBA multi-factor designs (multifactor skill).
- MOFA multi-view integration (multi-omics skill).
- LIPID MAPS lipid parsing + class aggregation (lipidomics skill).

If a future tutorial only swaps the imputation method, normalization choice, or adds a new transformation primitive, update this skill rather than creating a duplicate.

## Branch Selection

**Imputer (`impute(method=...)`) by missingness mechanism**
- `qrilc` — left-censored / below-detection (LC-MS default). Draws from a truncated quantile-regressed normal.
- `knn` — missing-at-random (sample mishandling, integration glitches). Sensitive to `n_neighbors`.
- `half_min` — fast baseline, fills with feature half-minimum. Use only when distributional fidelity is not critical.
- `zero` — only when zeros are structurally meaningful (e.g., truly absent metabolite); never the right call for LC-MS NaNs.

**Sample normalizer (`normalize(method=...)`)**
- `pqn` — Probabilistic Quotient Normalization. **Default for both NMR and LC-MS.**
- `tic` — Total Ion Count. Only safe when the matrix is conserved (all samples on the same instrument run, same injection volume).
- `median` — sample-median ratio. Robust to extreme features; less aggressive than PQN.
- `mstus` — MS Total Useful Signal (LC-MS only — drops low-coverage features before TIC).

**Feature transform (`transform(method=...)`) — sequence matters**
- `log` then downstream univariate stats (`differential` / `anova`). Pass `log_transformed=True` so fold-changes are computed on un-logged means.
- `log` → `pareto` for downstream multivariate (`plsda`, `opls_da`, `asca`). Pareto centres each feature and divides by `sqrt(SD)` — preserves dominant signal structure.
- `autoscale` (z-score) only if you specifically want unit variance per feature.
- `glog` for very-low-intensity LC-MS features where `log` over-amplifies noise.

**Correction stack (LC-MS only) — order matters**
1. `drift_correct` — injection-order drift via LOESS on QC. Skip if QC samples aren't interspersed.
2. `serrf` — within-batch QC variability via Random Forest. Needs ≥5 QC pools per batch.
3. `ov.bulk.batch_correction` (ComBat) — between-batch additive shifts. Use only if residual batch effect remains after SERRF and you have ≥3 samples per batch.

**`pyMetabo` vs. functional API**
- Use `pyMetabo` for reproducible single-seed pipelines and when downstream code wants `.deg_table` / `.plsda_result` / `.significant_metabolites()` accessors.
- Use the functional API when you want to inspect intermediate `AnnData` between every step or chain in / out of non-metabol code.

## Input Contract

- `AnnData` with `obs=samples`, `var=metabolites` (loaders enforce this — pass `transpose=True` if your file is the other way).
- `adata.obs` must carry the factor column you'll pass as `group_col`; the loader writes it as `obs['group']`.
- For LC-MS drift / SERRF: `adata.obs` needs `sample_type` (with a `'QC'` label) and `batch`; for `drift_correct` it also needs an injection-order column (numeric, monotonic).
- Missing values must be NaN, not zero (loaders handle this for the supported formats; check if you built the AnnData manually).
- `adata.X` is dense float — sparse matrices are not supported by the imputers / SERRF.

## Minimal Execution Patterns

```python
# Functional one-shot (NMR / MetaboAnalyst CSV)
import omicverse as ov

ov.plot_set()
adata = ov.metabol.read_metaboanalyst('human_cachexia.csv', group_col='Muscle loss')
adata = ov.metabol.impute(adata, method='qrilc', q=0.01, seed=0)
adata = ov.metabol.normalize(adata, method='pqn')
log_adata = ov.metabol.transform(adata, method='log')
pareto_adata = ov.metabol.transform(log_adata, method='pareto', stash_raw=False)

deg = ov.metabol.differential(
    log_adata, group_col='group', group_a='cachexic', group_b='control',
    method='welch_t', log_transformed=True,
)
fig, ax = ov.metabol.volcano(deg, padj_thresh=0.10, log2fc_thresh=0.3, label_top_n=8)
```

```python
# Chained class API — same numerics, more compact
m = (
    ov.metabol.pyMetabo(adata.copy())
      .impute(method='qrilc', seed=0)
      .normalize(method='pqn')
      .transform(method='log')
      .differential(group_col='group', group_a='cachexic', group_b='control',
                    method='welch_t', log_transformed=True)
      .transform(method='pareto', stash_raw=False)
)
m.deg_table.head()
m.significant_metabolites(padj_thresh=0.10, log2fc_thresh=0.3)
```

```python
# LC-MS correction stack
adata = ov.metabol.drift_correct(
    adata,
    injection_order='order',
    qc_mask=(adata.obs['sample_type'] == 'QC').to_numpy(),
    frac=0.5,
)
adata = ov.metabol.serrf(
    adata,
    qc_col='sample_type', qc_label='QC',
    batch_col='batch',
    top_k=8, n_estimators=100, seed=0,
)
ov.bulk.batch_correction(adata, batch_key='batch', key_added='batch_correction')
# Corrected matrix at adata.layers['batch_correction']

real = adata[adata.obs['sample_type'] == 'real'].copy()
qc_df = ov.metabol.sample_qc(real, n_components=3, alpha=0.95)
ov.metabol.sample_qc_plot(qc_df)
flagged = qc_df.index[qc_df['is_outlier']].tolist()
```

## Validation

- `'group' in adata.obs.columns` after loading — if missing, the loader's `group_col` was wrong.
- `not np.isnan(adata.X).any()` after `impute(...)` — if NaNs remain, missingness exceeded `missing_threshold` for at least one column; either lower the threshold or drop those features.
- After `normalize(method='pqn')`, row-median dispersion (`max/min`) should drop versus the raw matrix — if it doesn't, picked the wrong method.
- After `transform(method='pareto', ...)`, every feature's mean should be ≈0 — if not, Pareto wasn't applied (likely because data wasn't log-transformed first).
- After `serrf`, `adata.var['cv_qc_serrf']` should be lower than `cv_qc_raw` for most features; the scatter `cv_qc_raw` vs `cv_qc_serrf` should sit below the diagonal.
- After `sample_qc`, manually inspect samples flagged with `is_outlier=True` before dropping — biological outliers (rare phenotypes) and technical outliers look identical here.
- The `pyMetabo` chain and the functional API must agree to machine precision when the same seed and method choices are used (`assert (m.deg_table['pvalue'].values == deg['pvalue'].values).all()`).
- If you only validated the chain on the cachexia NMR dataset, do not claim LC-MS correctness; both ingest paths share imputation / normalization but diverge on `read_lcms` parsing and on the drift / SERRF / ComBat stack.

## Resource Map

- See [`reference.md`](reference.md) for copy-paste-ready snippets for each stage.
- See [`references/source-grounding.md`](references/source-grounding.md) for the live API signatures and the docstring-supplementation log (`pyMetabo` chainable methods + plotting helpers were docstring-empty; this skill drove their backfill).
- For multivariate (PLS-DA / OPLS-DA / VIP / S-plot) downstream of this preprocessing, see the **multivariate** skill (separate).
- For pathway enrichment (MSEA ORA / GSEA, mummichog), see the **pathway-multifactor** skill (separate).
- For biomarker selection with nested CV, see the **biomarker** skill.
- For multi-omics integration (MOFA on metabolomics + RNA-seq), see the **multi-omics MOFA** skill.

## Examples
- "Load `human_cachexia.csv` (MetaboAnalyst format, factor column `Muscle loss`), impute with QRILC, PQN-normalize, log-transform, run a Welch t-test for cachexic vs. control, and plot a volcano with padj<0.10."
- "My LC-MS run has injection-order drift and 3 batches with QC pools — apply drift_correct + SERRF + ComBat in that order; confirm SERRF dropped median QC CV% by at least 30 %."
- "Run sample QC on the real-sample subset and list any samples above the 0.95 critical Hotelling line."

## References
- Tutorial notebooks:
  - [`t_metabol_01_intro.ipynb`](https://omicverse.readthedocs.io/en/latest/Tutorials-metabol/t_metabol_01_intro/) — full pipeline on the cachexia NMR cohort.
  - [`t_metabol_06_batch_correction.ipynb`](https://omicverse.readthedocs.io/en/latest/Tutorials-metabol/t_metabol_06_batch_correction/) — drift / SERRF / ComBat on synthetic LC-MS.
- Live API verified against `omicverse.metabol` at the time of writing — see [`references/source-grounding.md`](references/source-grounding.md).
