# Source Grounding — Bulk Metabolomics Preprocessing

## Interfaces Checked

The skill content was grounded against live `omicverse.metabol` (v0.4.0) and `omicverse.bulk` source using:

- `inspect.signature(...)` on every public function listed in the Interface Summary.
- `inspect.getdoc(...)` on every public function plus the `pyMetabo` lifecycle class and its 11 chainable methods.
- Direct reading of `omicverse/metabol/{__init__,pymetabo,_qc,_batch,_impute,_norm,_transform,_stats,plotting}.py` to confirm parameter behavior, in-place mutations, and the layers / obs columns each function writes.
- Cross-checking the `t_metabol_01_intro.ipynb` and `t_metabol_06_batch_correction.ipynb` tutorial cells against those signatures — every API call cited in this skill appears in at least one tutorial cell or in the `omicverse.metabol` `__init__.py` `__all__` export list.

## Live signatures — loaders

```python
ov.metabol.read_metaboanalyst(
    path: 'str | Path', *,
    group_col: 'str',
    sample_col: 'Optional[str]' = None,
    transpose: 'bool' = False,
) -> 'AnnData'

ov.metabol.read_lcms(
    path: 'str | Path', *,
    feature_id_sep: 'str' = '/',
    sample_col: 'Optional[str]' = None,
    group_col: 'Optional[str]' = None,
    label_row: 'Optional[str]' = None,
    transpose: 'bool' = True,
) -> 'AnnData'

ov.metabol.read_wide(
    path: 'str | Path', *,
    sample_col: 'Optional[str]' = None,
    group_col: 'Optional[str]' = None,
    transpose: 'bool' = False,
) -> 'AnnData'
```

Source-grounded behavior:
- All loaders write `obs['group']` from the `group_col` argument so downstream stats see a uniform column name.
- `read_lcms` parses the feature-name token (split on `feature_id_sep`) into `var['mz']` and `var['rt']` numeric columns.
- `transpose=True` (LC-MS default) rotates rows-as-features files into rows-as-samples — the AnnData convention.

## Live signatures — preprocessing chain

```python
ov.metabol.impute(
    adata, *,
    method: ImputeMethod = 'qrilc',
    missing_threshold: float = 0.5,
    n_neighbors: int = 5,
    q: float = 0.01,
    seed: int = 0,
) -> AnnData

ov.metabol.normalize(
    adata, *,
    method: NormMethod = 'pqn',
    reference: Literal['median', 'mean'] = 'median',
    missing_threshold: float = 0.5,
) -> AnnData

ov.metabol.transform(
    adata, *,
    method: TransformMethod = 'log',
    pseudocount: float = 1.0,
    stash_raw: bool = True,
) -> AnnData
```

Source-grounded behavior:
- `impute` drops features whose missing fraction exceeds `missing_threshold` *before* imputation (otherwise QRILC's quantile estimate is unstable).
- `qrilc` uses a left-truncated normal at the `q`-th feature quantile to draw fills; `seed` gates reproducibility.
- `normalize(method='pqn')` computes the reference profile as the per-feature median (or mean) over all samples (post-row-sum normalization), then divides each sample by the median ratio.
- `transform(method='pareto')` does both centering and `1/sqrt(SD)` scaling — not just scaling. Confirmed in `_transform.py`.
- `stash_raw=True` writes the prior `adata.X` into `adata.layers['raw_X']`; the second call with `stash_raw=False` avoids overwriting it.

## Live signatures — QC / MS correction

```python
ov.metabol.cv_filter(
    adata, *,
    qc_mask: str | np.ndarray | None = None,
    cv_threshold: float = 0.30,
    across: str = 'qc',  # or 'all'
) -> AnnData

ov.metabol.blank_filter(
    adata, *,
    blank_mask: str | np.ndarray,
    ratio: float = 3.0,
) -> AnnData

ov.metabol.drift_correct(
    adata, *,
    injection_order: str | np.ndarray,
    qc_mask: str | np.ndarray,
    frac: float = 0.5,
) -> AnnData

ov.metabol.serrf(
    adata, *,
    qc_col: str,
    qc_label: str = 'QC',
    batch_col: Optional[str] = None,
    top_k: int = 10,
    n_estimators: int = 100,
    min_qc_samples: int = 5,
    layer: Optional[str] = None,
    seed: int = 0,
) -> AnnData

ov.metabol.sample_qc(
    adata, *,
    n_components: int = 2,
    alpha: float = 0.95,
    center: bool = True,
    scale: bool = True,
    layer: str | None = None,
) -> pd.DataFrame  # columns: T2, DModX, T2_crit, DModX_crit, is_outlier
```

Source-grounded behavior:
- `cv_filter(across='qc')` computes CV over QC samples only (the canonical SERRF preflight); `across='all'` over every sample.
- `drift_correct` uses `statsmodels.nonparametric.smoothers_lowess.lowess` per-feature; `frac` is forwarded directly.
- `serrf` adds two columns to `var`: `cv_qc_raw` (pre-correction) and `cv_qc_serrf` (post-correction). When `batch_col` is set, RFs are trained per batch independently — that's why batch must be present *before* SERRF, not after.
- `serrf` raises if any batch has fewer than `min_qc_samples=5` QC pools.
- `sample_qc` returns critical values from F-distribution (T²) and chi² (DModX) at level `alpha`; both critical values are written as constant columns to make plotting trivial.

ComBat (cross-module call):

```python
ov.bulk.batch_correction(
    adata, *,
    batch_key: str,
    key_added: str = 'batch_correction',
)
```

Source-grounded behavior:
- Writes the corrected matrix to `adata.layers[key_added]`; `adata.X` is left untouched.
- Implementation is the same `pyComBat`-style ComBat used by `bulk-combat-correction` skill; the metabol-side use case differs only in pre-conditioning (PQN + log).

## Live signatures — univariate stats + volcano

```python
ov.metabol.differential(
    adata, *,
    group_col: str = 'group',
    group_a: Optional[str] = None,
    group_b: Optional[str] = None,
    method: TestMethod = 'welch_t',
    layer: Optional[str] = None,
    log_transformed: bool = True,
) -> pd.DataFrame  # log2fc, pvalue, padj, mean_a, mean_b, ...

ov.metabol.anova(
    adata, *,
    group_col: str = 'group',
    groups: Optional[list] = None,
    method: ANOVAMethod = 'welch_anova',
    layer: Optional[str] = None,
) -> pd.DataFrame

ov.metabol.volcano(
    deg, *,
    padj_thresh: float = 0.05,
    log2fc_thresh: float = 1.0,
    label_top_n: int = 10,
    use_pvalue: bool = False,
    clip_log2fc: Optional[float] = None,
    ax: Optional[plt.Axes] = None,
    figsize: tuple[float, float] = (5.5, 4.5),
)
```

Source-grounded behavior:
- `log_transformed=True` instructs `differential` to undo the log before computing fold-changes, so `log2fc` is in the canonical sense (positive ⇒ higher in group_a). Skipping this flag on log-data inflates `log2fc` magnitudes.
- BH-FDR is the only multiple-testing correction; for multi-test family planning use `padj` directly.

## Live signatures — `pyMetabo` chainable class

```python
class pyMetabo:
    adata: AnnData
    random_state: int = 0
    raw: AnnData          # frozen copy of input (provenance)
    deg_table: Optional[pd.DataFrame]
    plsda_result: Optional[PLSDAResult]

    # all return self; all keyword-only
    .cv_filter(qc_mask, cv_threshold=0.30)
    .drift_correct(injection_order, qc_mask, frac=0.5)
    .blank_filter(blank_mask, ratio=3.0)
    .impute(method='qrilc', seed=None, **kwargs)            # falls back to self.random_state
    .normalize(method='pqn', **kwargs)
    .transform(method='log', **kwargs)
    .differential(group_col='group', group_a, group_b, method='welch_t', log_transformed=True)
    .plsda(n_components=2, group_col='group', ...)
    .opls_da(n_ortho=1, group_col='group', ...)
    # consumers
    .vip_table() -> pd.DataFrame
    .significant_metabolites(padj_thresh=0.05, log2fc_thresh=1.0) -> pd.DataFrame
```

## Docstring supplementation log

The following symbols had **empty or 1-line** docstrings as observed during this skill's source-grounding pass. They have been filled with proper Numpy-style documentation (parameters, return values, biological/statistical rationale, ordering notes) as part of producing this skill. Diff applies to `omicverse/metabol/pymetabo.py` and `omicverse/metabol/plotting.py`.

| Symbol | Prior state | Action |
|---|---|---|
| `pyMetabo.cv_filter` | empty | filled — points at functional `cv_filter`, calls out `across='qc'` semantics |
| `pyMetabo.drift_correct` | empty | filled — explains LOESS-on-QC and injection-order requirement |
| `pyMetabo.blank_filter` | empty | filled — explains LC-MS blank-contaminant pruning |
| `pyMetabo.impute` | empty | filled — guidance on `qrilc` (MNAR) vs `knn` (MAR) vs `half_min` / `zero` baselines; reproducibility behavior |
| `pyMetabo.normalize` | empty | filled — PQN canonical; ordering relative to imputation/transform |
| `pyMetabo.transform` | empty | filled — `log → pareto` ordering for univariate→multivariate |
| `pyMetabo.differential` | empty | filled — `log_transformed` semantics; result schema |
| `pyMetabo.plsda` | empty | filled — Pareto-prerequisite; downstream `vip_table()` |
| `pyMetabo.opls_da` | empty | filled — orthogonal-component motivation; cleaner S-plot vs PLS-DA |
| `pyMetabo.vip_table` | empty | filled — VIP > 1 threshold; precondition |
| `pyMetabo.significant_metabolites` | empty | filled — relaxation guidance for small-cohort metabolomics |
| `metabol.plotting.vip_bar` | 1-line | expanded — colour-by-coef-sign explained, VIP=1 cutoff, parameter table |
| `metabol.plotting.sample_qc_plot` | 1-line | expanded — SIMCA-style T²-vs-DModX semantics |
| `metabol.plotting.dgca_class_bar` | 1-line | expanded — DC-class encoding, palette rationale |

Other public symbols (`read_metaboanalyst`, `impute`, `normalize`, `transform`, `differential`, `anova`, `drift_correct`, `serrf`, `sample_qc`, `volcano`, `s_plot`, `pathway_bar`, `pathway_dot`, `asca`, `mixed_model`, `meba`, `roc_feature`, `biomarker_panel`, `dgca`, `run_mofa`, `msea_ora`, `msea_gsea`, `mummichog_basic`, `annotate_peaks`, `parse_lipid`, `annotate_lipids`, `aggregate_by_class`, `lion_enrichment`, `cv_filter`, `blank_filter`) already had ≥10-line Numpy-style docstrings in the package — verified, no changes required.

## Reviewer-Run Empirical Checks

- `omicverse` v from `inspect.getmodule` matches the local checkout at `/scratch/users/steorra/analysis/omicverse_dev/omicverse/` (commit `78586cbe`+ on `master`).
- Notebook calls in `t_metabol_01_intro.ipynb` execute the documented signatures as-is (verified by reading the code cells, not by running them).
- Notebook calls in `t_metabol_06_batch_correction.ipynb` use synthetic data construction inline; the LC-MS-correction signatures match the source.
- No live smoke run was executed for this skill (per the project's "skip smoke when no validated dataset is on hand" policy).
