# Source Grounding — Paired microbiome × metabolomics

## Interfaces Checked

`omicverse.micro._pair` — `paired_spearman`, `paired_cca`, `MMvec`, `simulate_paired`, `fetch_franzosa_ibd_2019`, plus the five plotting helpers (`plot_mmvec_training`, `plot_cca_scatter`, `plot_cooccurrence`, `plot_embedding_biplot`, `plot_paired_method_comparison`). Verified via `inspect.signature` + `inspect.getdoc` + direct reading of `omicverse/micro/_pair.py`. Cell-by-cell cross-checking against `t_micro_metabol_paired.ipynb`.

## Live signatures

```python
ov.micro.fetch_franzosa_ibd_2019(
    data_dir: str,
    overwrite: bool = False,
    microbe_count_scale: float = 1e6,
) -> Tuple[ad.AnnData, ad.AnnData]

ov.micro.simulate_paired(
    n_samples: int = 30, n_microbes: int = 40,
    n_metabolites: int = 20, n_pairs: int = 5,
    effect_range: Tuple[float, float] = (1.0, 2.0),
    depth_range: Tuple[int, int] = (1000, 10000),
    seed: int = 0,
) -> Tuple[ad.AnnData, ad.AnnData, pd.DataFrame]

ov.micro.paired_spearman(
    adata_microbe: ad.AnnData,
    adata_metabolite: ad.AnnData,
    clr_microbe: bool = True,
    log1p_metabolite: bool = True,
    min_prevalence: float = 0.0,
) -> pd.DataFrame

ov.micro.paired_cca(
    adata_microbe: ad.AnnData,
    adata_metabolite: ad.AnnData,
    n_components: int = 3,
    clr_microbe: bool = True,
    log1p_metabolite: bool = True,
    max_iter: int = 500,
) -> Dict[str, Any]

class MMvec:
    def __init__(self, n_latent: int = 3, lr: float = 0.05,
                 epochs: int = 1000, val_frac: float = 0.1,
                 patience: int = 100, l2: float = 1e-3,
                 seed: int = 0, device: Optional[str] = None): ...
    def fit(self, adata_microbe, adata_metabolite, verbose=False) -> 'MMvec': ...
    def cooccurrence(self) -> pd.DataFrame: ...
    def conditional_probabilities(self) -> pd.DataFrame: ...
    def top_pairs(self, n: int = 20) -> pd.DataFrame: ...
    @property
    def microbe_embeddings_(self) -> pd.DataFrame: ...
    @property
    def metabolite_embeddings_(self) -> pd.DataFrame: ...
    # populated by fit:
    #   U_, V_, beta_, microbe_names_, metabolite_names_,
    #   loss_history_, val_loss_history_, best_epoch_

ov.micro.plot_mmvec_training(mmvec, ax=None)
ov.micro.plot_cca_scatter(cca_result, component=1, ax=None)
ov.micro.plot_cooccurrence(score_df, top_n=15, ax=None, cmap='RdBu_r')
ov.micro.plot_embedding_biplot(mmvec, components=(0, 1), label_top=5, ax=None)
ov.micro.plot_paired_method_comparison(truth, spearman_df=None, mmvec_model=None, ax=None)
```

## Source-grounded behavior

**`fetch_franzosa_ibd_2019`:**
- Downloads + parses two CSVs from the original Franzosa 2019 paired IBD dataset; caches under `data_dir`. Both AnnDatas share `obs_names` and have `obs['Study.Group']` (Crohn / UC / Control).
- Microbe counts are stored as integers (after a per-microbe `microbe_count_scale=1e6` × proportion conversion); metabolites are intensities. `adata_mt.var['name']` carries human-readable metabolite names where available.

**`paired_spearman`:**
- Internally computes Spearman ρ per (microbe, metabolite) pair using `scipy.stats.spearmanr`. Pre-applies CLR to microbes (`clr_microbe=True`) and `log1p` to metabolites (`log1p_metabolite=True`) by default.
- `min_prevalence` filters microbes (and only microbes) below the threshold before testing. Returns one row per surviving (microbe, metabolite) pair with `r`, `pvalue`, `fdr_bh`.

**`paired_cca`:**
- Wraps `sklearn.cross_decomposition.CCA(n_components, max_iter)`. Returns a dict (not a fitted object) with `canonical_correlations` (length `n_components`), `microbe_scores` / `metabolite_scores` (samples × components), `microbe_loadings` / `metabolite_loadings` (features × components).
- `_check_paired` enforces `obs_names` alignment.

**`MMvec.fit`:**
- Implementation is ~80 lines of PyTorch. Loss = per-microbe negative log-likelihood of the sample-weighted metabolite distribution.
- Cooccurrence weights `W[i, j] = Σ_s X_mb[s, i] · X_mt[s, j] / Σ_s X_mt[s, :]`.
- Logits `U @ V.T + β` softmaxed per row; loss is `-(W * log_softmax(logits)).sum() / W.sum()`.
- Train/val split by `val_frac` of samples; cooccurrence weights computed separately on each split; early-stop on the validation set with `patience` epochs of no improvement.
- Best-checkpoint state is restored at end of training; if no validation set (val_frac=0), final-epoch state is kept.
- `device='cuda'` auto-selected when available; falls back to CPU otherwise.

**`MMvec.cooccurrence`:** symmetric `U @ V.T` (no softmax). Use for top-pair queries.

**`MMvec.conditional_probabilities`:** row-softmax of `U @ V.T + β` — proper P(metabolite | microbe). Numerically stable (subtracts row max).

**`MMvec.top_pairs(n)`:** stacks `cooccurrence()` to long format, sorts by `|score|`, returns top `n` with columns `microbe`, `metabolite`, `score`. Score sign is preserved so caller can see direction.

**`simulate_paired`:**
- Plants `n_pairs` (microbe, metabolite) pairs with effect drawn from `effect_range`. Samples receive variable depth from `depth_range`. Returns (`adata_mb`, `adata_mt`, `truth_df`) where `truth_df` has columns `microbe`, `metabolite` for the planted pairs.

**Plotting helpers** — all backed by source-grounded reading + cross-checked with notebook usage.

## Notebook ↔ skill alignment

| Notebook section | Skill section |
|---|---|
| `fetch_franzosa_ibd_2019` ingest + group breakdown | Quick Workflow §1 |
| Top-150 microbes / top-200 annotated metabolites filter | Quick Workflow §2 |
| `paired_spearman(min_prevalence=0.1)` + FDR count | Quick Workflow §3 |
| `paired_cca(n_components=3)` + canonical correlations | Quick Workflow §4 |
| `plot_cca_scatter(component=1)` | Quick Workflow §5 |
| `MMvec(n_latent=4, epochs=600, val_frac=0.15, patience=80, seed=0).fit(...)` + best_epoch / loss | Quick Workflow §6 |
| `plot_mmvec_training`, `plot_cooccurrence`, `plot_embedding_biplot` | Quick Workflow §7 |
| Synthetic recovery: `simulate_paired(n_pairs=5)` + Spearman + MMvec + `plot_paired_method_comparison` | Quick Workflow §8; reference.md synthetic block |

## Docstring supplementation log

Driven by this skill's source-grounding pass. Diff applies to `omicverse/micro/_pair.py`:

| Symbol | Prior state | Action |
|---|---|---|
| `MMvec.fit` | 1-line | filled — full Parameters + algorithm sketch (cooccurrence weights, logits formula, loss, early-stop / best-checkpoint behavior). |
| `MMvec.cooccurrence` | 1-line | filled — explains symmetric vs row-normalised distinction; pairs with `top_pairs` and `conditional_probabilities`. |
| `MMvec.conditional_probabilities` | 1-line | filled — explains P(metabolite \| microbe) semantics + numerical-stability detail. |
| `MMvec.top_pairs` | 1-line | filled — explains absolute-score sort and sign-preservation. |
| `plot_mmvec_training` | 1-line | filled — describes loss curve interpretation (overfit vs healthy). |
| `plot_embedding_biplot` | 1-line | filled — describes biplot semantics + `label_top` parameter. |

Other functions in `_pair.py` had complete docstrings (`fetch_franzosa_ibd_2019` 37L, `paired_spearman` 18L, `paired_cca` 10L, the `MMvec` class itself 34L, `simulate_paired` 5L which is short but adequate, `plot_cca_scatter` 10L, `plot_cooccurrence` 5L, `plot_paired_method_comparison` 4L). Left as-is.

## Reviewer-Run Empirical Checks

- All cited functions importable: `from omicverse.micro import (paired_spearman, paired_cca, MMvec, simulate_paired, fetch_franzosa_ibd_2019, plot_mmvec_training, plot_cca_scatter, plot_cooccurrence, plot_embedding_biplot, plot_paired_method_comparison)` ✓
- `MMvec.fit` early-stop algorithm: validates that `best_epoch_` is set when val_frac > 0 and `patience_left` decrements correctly (read from source).
- Tutorial cells use exactly the verified signatures.
- No live smoke run executed — Franzosa download + MMvec training (~5 min) is the canonical smoke target; `simulate_paired(n_pairs=5)` + MMvec + Spearman is a second-tier acceptance probe.
