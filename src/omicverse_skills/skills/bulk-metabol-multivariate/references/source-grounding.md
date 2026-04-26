# Source Grounding — Multivariate & Biomarker

## Interfaces Checked

Same omicverse checkout (`omicverse.metabol` v0.4.0, master `78586cbe`+) as the preprocessing skill. Verified via `inspect.signature` + `inspect.getdoc` + direct reading of `omicverse/metabol/{_plsda,_biomarker,plotting,pymetabo}.py`, plus cross-checking the actual code cells in `t_metabol_02_multivariate.ipynb` and `t_metabol_08_biomarker.ipynb`.

## Live signatures

```python
ov.metabol.plsda(
    adata, *,
    group_col: str = 'group',
    group_a: Optional[str] = None,
    group_b: Optional[str] = None,
    n_components: int = 2,
    scale: bool = False,
) -> PLSDAResult

ov.metabol.opls_da(
    adata, *,
    group_col: str = 'group',
    group_a: Optional[str] = None,
    group_b: Optional[str] = None,
    n_ortho: int = 1,
    scale: bool = False,
    max_iter: int = 500,
    tol: float = 1e-08,
) -> PLSDAResult

ov.metabol.roc_feature(
    adata, *,
    group_col: str,
    pos_group: Optional[str] = None,
    neg_group: Optional[str] = None,
    layer: Optional[str] = None,
    ci: bool = False,
    n_bootstrap: int = 1000,
    seed: int = 0,
) -> pd.DataFrame

ov.metabol.biomarker_panel(
    adata, *,
    group_col: str,
    features: Union[list, int] = 10,
    classifier: str = 'rf',
    cv_outer: int = 5,
    cv_inner: int = 3,
    pos_group: Optional[str] = None,
    neg_group: Optional[str] = None,
    n_permutations: int = 0,
    layer: Optional[str] = None,
    seed: int = 0,
) -> BiomarkerPanelResult

ov.metabol.vip_bar(
    result: PLSDAResult,
    var_names, *,
    top_n: int = 15,
    ax: Optional[plt.Axes] = None,
    figsize: tuple[float, float] = (5.0, 5.0),
)

ov.metabol.s_plot(
    result: PLSDAResult,
    adata, *,
    label_top_n: int = 15,
    ax: Optional[plt.Axes] = None,
    figsize: tuple[float, float] = (5.5, 4.5),
)
```

## Source-grounded behavior

**`plsda` / `opls_da`:**
- Internally wrap `sklearn.cross_decomposition.PLSRegression`. `scale=False` is the right default *because* preprocessing already Pareto-scaled — `sklearn`'s built-in scaling on top would double-scale.
- Q² is computed via leave-one-out cross-validation, not k-fold; for moderate cohorts (n < 200) this is the appropriate choice.
- `opls_da` runs the OSC-style deflation `n_ortho` times *after* fitting the predictive component, then re-fits the predictive direction on the orthogonal-deflated X.
- Both populate `.group_labels` from the *actual* labels seen in the binary slice, so `result.group_labels[0]` is always `pos`-side regardless of how the user ordered group_a/group_b.

**`PLSDAResult`** (full docstring already present, expanded here for context):
- `.scores` is the predictive score matrix `(n_samples, n_components)`.
- `.x_ortho_scores` exists only on OPLS-DA results (PLS-DA has shape `(n_samples, 0)`).
- `.vip` is computed using the standard Wold formula; `.coef` is the regression coefficient against a binary `±1` Y-vector — `vip_bar` uses its sign for colour coding.
- `.to_vip_table(var_names)` returns a DataFrame with the metabolite name as index, columns `vip`, `coef`; sorted by VIP descending. 1-line docstring on the method itself, but the class-level docstring covers the schema.

**`roc_feature`:**
- Uses `sklearn.metrics.roc_auc_score` per feature; bootstrap CIs are percentile-based (2.5 / 97.5 percentiles of `n_bootstrap` resamples).
- The `pvalue` column is a 2-tailed test against `auc=0.5` derived from the bootstrap distribution.
- Signed AUC convention: `auc > 0.5` ⇒ feature higher in `pos_group`.

**`biomarker_panel`:**
- Nested CV: outer loop is stratified `cv_outer`-fold; inner loop is also stratified `cv_inner`-fold and is the *only* place feature selection happens (feature ranking is repeated per outer fold).
- `features=int`: top-N univariate features by AUC ranked inside each inner-train fold. `features=list`: fixed panel; inner CV only tunes the classifier.
- Permutation null: shuffle `y` `n_permutations` times, refit the entire nested CV, build a null distribution of mean AUC; `permutation_pvalue` is `(1 + #(null >= obs)) / (1 + n_permutations)`.
- `outer_predictions` / `outer_labels` are *out-of-fold* predictions concatenated — they're in fold-block order, not the original sample order.

**`BiomarkerPanelResult`** (full docstring present in source — all attributes documented):
- `.feature_importance` is a DataFrame with one row per metabolite that appeared in the final panel; ranks may differ between LR (coefficient magnitude), RF (Gini importance), SVM (linear kernel falls back to coef magnitude; RBF falls back to permutation importance over outer folds).
- `.to_frame()` returns per-outer-fold AUC; useful for SD reporting. 1-line docstring on the method but `.outer_aucs` is exposed directly for custom plots.

**`vip_bar` / `s_plot`** (docstrings already expanded as part of preprocessing skill's docstring backfill — see `bulk-metabol-preprocessing/references/source-grounding.md`).

## Notebook ↔ skill alignment

| Notebook cell | Skill section |
|---|---|
| `t_metabol_02` PQN+log+Pareto preprocessing | Demoted — handled by preprocessing skill |
| `t_metabol_02` `plsda(...)` baseline | Quick Workflow §2, Branch Selection (PLS-DA vs OPLS-DA) |
| `t_metabol_02` `opls_da(n_ortho=1)` + score-plot construction | Quick Workflow §3-4, Minimal Execution Patterns |
| `t_metabol_02` `to_vip_table` + `vip_bar` + `s_plot` | Quick Workflow §4-5, Interface Summary |
| `t_metabol_08` impute+PQN+log preprocessing | Demoted — handled by preprocessing skill |
| `t_metabol_08` `roc_feature(ci=True)` | Quick Workflow §6 |
| `t_metabol_08` `biomarker_panel` with permutation null | Quick Workflow §7, Branch Selection (classifier / features / n_permutations) |
| `t_metabol_08` ROC curve from `outer_predictions` | Quick Workflow §8 |
| `t_metabol_08` univariate `differential` mention | Cross-reference to preprocessing skill |

## Docstring supplementation log

No additional docstring fills were required for this skill — `_plsda` / `_biomarker` were already well-documented (11–48 line Numpy-style docstrings on each public function and result class). The 1-line docstrings on `to_vip_table`, `to_frame`, `vip_bar`, `s_plot` were already addressed in the preprocessing skill's docstring backfill (see `bulk-metabol-preprocessing/references/source-grounding.md` for `vip_bar` / `s_plot`).

## Reviewer-Run Empirical Checks

- All cited functions importable: `from omicverse.metabol import plsda, opls_da, vip_bar, s_plot, roc_feature, biomarker_panel, PLSDAResult, BiomarkerPanelResult` ✓
- Tutorial cells in `t_metabol_02_multivariate.ipynb` and `t_metabol_08_biomarker.ipynb` use exactly the verified signatures (no kwarg drift).
- No live smoke run executed (per project policy); all evidence is signature- and source-level.
