---
name: omicverse-bulk-celltype-deconvolution
description: Cell-type composition deconvolution of bulk RNA-seq using a single-cell reference. Wraps `ov.bulk.Deconvolution` with four interchangeable backends — TAPE, Scaden, BayesPrism, OmicsTweezer — under one constructor / one `.deconvolution(method=...)` call. Use when inferring cell-type fractions from bulk samples given a paired scRNA-seq atlas, when reproducing `t_decov_bulk`, or when comparing deep-learning vs. Bayesian deconvolution methods on the same cohort.
---

# OmicVerse Bulk RNA-seq — Cell-Type Deconvolution

## Goal

Take a bulk RNA-seq cohort (`AnnData`, samples × genes) plus a single-cell reference (`AnnData`, cells × genes, with cell-type labels) and infer **cell-type fractions per bulk sample**. One unified class wraps four established methods under a `method=...` switch:

- **TAPE** (Cao 2022) — deep-learning autoencoder; default; CPU-friendly.
- **Scaden** (Menden 2020) — deep neural network; trains on pseudo-bulk mixtures; benefits from GPU.
- **BayesPrism** (Chu 2022) — full Bayesian model; multi-core CPU; produces posterior fractions.
- **OmicsTweezer** — joint reference-correction deconvolution.

Output is a `pd.DataFrame` indexed by bulk sample with one column per cell type — directly stackable as a per-sample bar chart, or summarisable per phenotype group via `ov.pl.plot_grouped_fractions`.

This skill is the **opposite direction** of `bulk-to-single-deconvolution` (Bulk2Single, which generates synthetic single cells from bulk). They share the word "deconvolution" but solve inverse problems — don't confuse them.

## Quick Workflow

1. Load bulk expression + single-cell reference, both as `AnnData`. The COVID-19 PBMC demo ships with both: `bulk_ad = ov.datasets.decov_bulk_covid_bulk()`, `single_ad_ref = ov.datasets.decov_bulk_covid_single()`.
2. Eyeball the reference: UMAP coloured by `cell.type.coarse` / `cell.type.fine` and by phenotype (`disease`); per-group composition with `ov.pl.cellproportion(adata=single_ad_ref, celltype_clusters='cell.type.coarse', groupby='disease')`.
3. **Construct**: `deconv_obj = ov.bulk.Deconvolution(adata_bulk=bulk_ad, adata_single=single_ad_ref, max_single_cells=10000, celltype_key='cell.type.coarse', cellstate_key='cell.type.fine', gpu=0)`. The `cellstate_key` (sub-celltype) is needed by some backends (BayesPrism's hierarchical model) and ignored by others.
4. **Run BayesPrism** (multi-core CPU, full Bayesian): `res = deconv_obj.deconvolution(method='bayesprism', n_cores=8, fast_mode=True)` → `pd.DataFrame` (samples × cell types).
5. **Run Scaden** (GPU-friendly DNN): re-instantiate with `gpu='mps'` (Apple) / `gpu=0` (CUDA index): `res2 = deconv_obj.deconvolution(method='scaden', scaler='ss', scale=True, datatype='counts', pseudobulk_size=2000)`.
6. **Re-order columns** to the canonical celltype order from the reference: `res = res[single_ad_ref.obs['cell.type.coarse'].cat.categories]`.
7. **Plot stacked bars**: build a colour dict from the reference's `uns['cell.type.coarse_colors']`, then `res.plot(kind='bar', stacked=True, figsize=(12, 3), color=color_dict)`.
8. (Optional) **Group-aware summary** with `ov.pl.plot_grouped_fractions(res, obs=bulk_ad.obs, group_key='disease', color_dict=color_dict)` — averages fractions per phenotype group with normalisation.
9. Compare two methods on the same cohort: build both `res` (BayesPrism) and `res2` (Scaden) DataFrames, plot side-by-side; the per-sample correlation reveals method-specific biases.

## Interface Summary

```python
ov.bulk.Deconvolution(
    adata_bulk: AnnData,
    adata_single: AnnData,
    max_single_cells: int = 5000,
    celltype_key: str = 'celltype',
    cellstate_key: str | None = None,
    gpu: int | str = 0,
)
# Constructor: subsamples adata_single if larger than max_single_cells,
# builds the signature matrix from celltype_key (and cellstate_key for
# hierarchical methods).

deconv_obj.deconvolution(
    method: str = 'tape',          # 'tape' | 'scaden' | 'bayesprism' | 'omicstweezer'
    sep: str = '\t',
    scaler: str = 'mms',           # TAPE / Scaden input scaling
    datatype: str = 'counts',
    genelenfile: str | None = None,
    mode: str = 'overall',         # TAPE running mode
    adaptive: bool = True,         # TAPE adaptive feature selection
    variance_threshold: float = 0.98,
    save_model_name: str | None = None,
    batch_size: int = 128,         # NN-based methods
    epochs: int = 128,             # NN-based methods
    seed: int = 1,
    scale_size: int = 2,
    scale: bool = True,
    n_cores: int = 4,              # BayesPrism only
    fast_mode: bool = True,        # BayesPrism approx updates
    pseudobulk_size: int = 2000,   # Scaden / TAPE training mixtures
    **kwargs,
) -> pd.DataFrame                  # samples × cell types
```

Built-in datasets:
- `ov.datasets.decov_bulk_covid_bulk() → AnnData` — Decov *et al.* 2020 COVID-19 PBMC bulk RNA-seq.
- `ov.datasets.decov_bulk_covid_single() → AnnData` — paired single-cell reference with `cell.type.coarse`, `cell.type.fine`, `disease` annotations.

Plotting helpers:
- `ov.pl.cellproportion(adata, celltype_clusters, groupby, legend=True, ax=None)` — stacked-bar of celltype proportions per group; reused for both reference QC (cells / disease) and downstream presentation of inferred fractions.
- `ov.pl.plot_grouped_fractions(res, obs, group_key, color_dict=None, agg='mean', normalize=True, figsize=(4, 4))` — phenotype-group-averaged stacked bars.

## Boundary

**Inside scope:**
- All four `ov.bulk.Deconvolution` backends (TAPE, Scaden, BayesPrism, OmicsTweezer).
- Reference subsampling via `max_single_cells`.
- Hierarchical celltype + cellstate handling (BayesPrism).
- Stacked-bar visualisation per sample / per phenotype group.
- Method comparison on the same cohort.

**Outside scope — separate skill:**
- Bulk2Single (generating synthetic single cells from bulk; the opposite direction) — see `bulk-to-single-deconvolution`.
- Spatial transcriptomics deconvolution (Tangram / cell2location / starfysh / FlashDeconv) — see `spatial-tutorials`.
- Building the upstream single-cell reference (preprocessing, annotation) — see existing single-cell skills.
- TCGA bulk preprocessing — see `tcga-preprocessing`.
- Differential expression on the bulk side — see `bulk-deg-analysis` / `bulk-deseq2-analysis`.

## Branch Selection

**Method by compute environment + interpretation needs**

| method        | Engine    | Compute       | Strengths                                   | Weaknesses                              |
|---------------|-----------|---------------|---------------------------------------------|-----------------------------------------|
| `tape`        | DNN (AE)  | CPU OK        | Default; sensible on most cohorts; adaptive feature selection | Less transparent than Bayesian methods  |
| `scaden`      | DNN       | GPU helps     | Robust on diverse cohorts; large benchmark history | Slower on CPU; trains per-cohort        |
| `bayesprism`  | Bayesian  | Multi-core CPU | Posterior fractions; hierarchical; handles low-quality references | No GPU acceleration; slowest on large cohorts |
| `omicstweezer`| Joint     | CPU OK        | Reference-correction integrated             | Newer; less benchmark history           |

**Default rule of thumb**: report at least **two** methods (e.g. BayesPrism + Scaden) and inspect agreement. Per-sample correlation between methods >0.8 across cell types is a healthy floor; <0.5 indicates the deconvolution is unstable for that sample (often due to a cell type that's poorly represented in the reference).

**`max_single_cells`**
- 5000 (default) — fast; appropriate for a 50k-cell reference (10 % subsample).
- 10000 (tutorial) — better resolution at moderate compute cost.
- Above 50000 — diminishing returns; BayesPrism in particular slows quadratically.

**`celltype_key` vs `cellstate_key`**
- `celltype_key`: required (the coarse annotation; output column granularity).
- `cellstate_key`: optional; when set, BayesPrism builds a hierarchical signature with sub-states inside each celltype. Improves accuracy when within-celltype heterogeneity is biologically meaningful.
- For TAPE / Scaden / OmicsTweezer, `cellstate_key` is ignored; only `celltype_key` matters.

**`gpu`**
- `gpu=0` (CUDA device 0) — default on a CUDA box; needed for fast Scaden / TAPE.
- `gpu='mps'` — Apple-Silicon GPU; works for Scaden in the tutorial.
- `gpu=-1` — force CPU; only practical for BayesPrism (which is CPU-native).
- BayesPrism ignores this — always CPU; tune `n_cores` instead.

**BayesPrism-specific knobs**
- `n_cores=8` — parallel chains; scale to physical cores.
- `fast_mode=True` — approximate updates; ~10x faster than full Gibbs with minor accuracy loss. Use `fast_mode=False` only for publication-grade inference when runtime allows.

**Scaden / TAPE-specific knobs**
- `scaler`: `'ss'` (StandardScaler) vs `'mms'` (MinMaxScaler — TAPE default). For Scaden the tutorial uses `'ss'`; for TAPE keep `'mms'`.
- `pseudobulk_size=2000` — number of synthetic mixtures generated to train the DNN. 2000 is fine for 5–15 cell types; raise to 5000+ for 30+ cell types.
- `batch_size=128`, `epochs=128` — defaults work on cohorts up to a few hundred bulk samples; raise epochs for larger cohorts.
- `datatype='counts'` (Scaden expects raw counts) vs `'tpm'` (TAPE accepts either).

**Re-ordering output columns** is a recurring tutorial pattern:
```python
res = res[single_ad_ref.obs[celltype_key].cat.categories]
```
Cell types in the result come back in arbitrary order; reordering matches the reference's canonical category order so the colour dict and stacked-bar legend stay consistent.

## Input Contract

- `adata_bulk`: `AnnData` with samples in `obs`, genes in `var`. Counts or normalised expression both work; the wrapper handles scaling per-method.
- `adata_single`: `AnnData` with `obs[celltype_key]` populated and (optionally) `obs[cellstate_key]`. Gene names in `var_names` must overlap with `adata_bulk.var_names` — the wrapper uses the intersection.
- For `BayesPrism`: install with `pip install BayesPrism` (or via the bundled extras).
- For `Scaden` / `TAPE`: bundled in OmicVerse via `omicverse.external`.
- `cell.type.coarse_colors` (or whatever the celltype palette key is) in `single_ad_ref.uns` — used to build the stacked-bar colour dict. Run `sc.pl.umap(adata, color=celltype_key)` once with scanpy if the colours aren't already populated; that writes the `_colors` array.

## Minimal Execution Patterns

```python
import omicverse as ov
ov.plot_set()

# 1) Built-in COVID-19 PBMC demo (replace with your own AnnData)
bulk_ad        = ov.datasets.decov_bulk_covid_bulk()
single_ad_ref  = ov.datasets.decov_bulk_covid_single()

# 2) Reference QC visuals
ov.pl.embedding(single_ad_ref, basis='X_umap',
                color=['cell.type.fine', 'cell.type.coarse'], ncols=1)
ov.pl.embedding(single_ad_ref, basis='X_umap', color='disease')

fig, ax = ov.plt.subplots(figsize=(1.5, 3))
ov.pl.cellproportion(
    adata=single_ad_ref,
    celltype_clusters='cell.type.coarse',
    groupby='disease', legend=True, ax=ax,
)
```

```python
# 3) BayesPrism — multi-core CPU, full Bayesian
deconv_obj = ov.bulk.Deconvolution(
    adata_bulk=bulk_ad,
    adata_single=single_ad_ref,
    max_single_cells=10000,
    celltype_key='cell.type.coarse',
    cellstate_key='cell.type.fine',
)
res = deconv_obj.deconvolution(
    method='bayesprism',
    n_cores=8, fast_mode=True,
)
res = res[single_ad_ref.obs['cell.type.coarse'].cat.categories]   # canonical order
print(res.head())                                                  # samples × cell types

# Stacked-bar visualisation
color_dict = dict(zip(
    single_ad_ref.obs['cell.type.coarse'].cat.categories,
    single_ad_ref.uns['cell.type.coarse_colors'],
))
ax = res.plot(kind='bar', stacked=True, figsize=(12, 3), color=color_dict)
ax.set(xlabel='Sample', ylabel='Cell Fraction',
       title='BayesPrism predicted fractions')
ov.plt.legend(bbox_to_anchor=(1.05, 1), ncol=1)
ov.plt.show()
```

```python
# 4) Scaden — GPU-friendly DNN; re-instantiate with gpu='mps' / 0 / -1
deconv_obj = ov.bulk.Deconvolution(
    adata_bulk=bulk_ad,
    adata_single=single_ad_ref,
    max_single_cells=10000,
    celltype_key='cell.type.coarse',
    cellstate_key='cell.type.fine',
    gpu='mps',                                # or gpu=0 for CUDA 0
)
res2 = deconv_obj.deconvolution(
    method='scaden',
    scaler='ss', scale=True, datatype='counts',
    pseudobulk_size=2000,
)
res2 = res2[single_ad_ref.obs['cell.type.coarse'].cat.categories]

ax = res2.plot(kind='bar', stacked=True, figsize=(12, 3), color=color_dict)
ax.set(xlabel='Sample', ylabel='Cell Fraction',
       title='Scaden predicted fractions')
ov.plt.legend(bbox_to_anchor=(1.05, 1), ncol=1)
ov.plt.show()
```

```python
# 5) Group-averaged comparison (e.g. healthy vs disease)
ov.pl.plot_grouped_fractions(
    res,
    obs=bulk_ad.obs,
    group_key='disease',
    color_dict=color_dict,
    agg='mean', normalize=True,
    figsize=(4, 4),
)
```

```python
# 6) Cross-method agreement check
import numpy as np
common = res.columns.intersection(res2.columns)
per_sample_corr = np.array([
    np.corrcoef(res.loc[s, common].values,
                res2.loc[s, common].values)[0, 1]
    for s in res.index
])
print(f'BayesPrism vs Scaden per-sample correlation: '
      f'{per_sample_corr.mean():.3f} ± {per_sample_corr.std():.3f}')
print(f'Samples with corr < 0.5: {(per_sample_corr < 0.5).sum()} / {len(per_sample_corr)}')
```

## Validation

- After construction: `deconv_obj.adata_single.n_obs <= max_single_cells`. The wrapper subsamples; if it doesn't shrink, the reference was already small enough.
- After `.deconvolution(...)`: result rows sum to ~1 (approximate compositional constraint). Violations of more than a few % are a red flag — re-check `celltype_key` mapping.
- After re-ordering columns: `set(res.columns) == set(single_ad_ref.obs[celltype_key].cat.categories)`. If columns are missing, the reference had a category with zero cells.
- Method agreement: per-sample correlation between two methods should mostly be >0.8 across cell types. Per-cell-type bias (e.g. BayesPrism systematically reports more T cells than Scaden) is a known method-level artefact; report it explicitly rather than picking the more "favourable" method.
- Cohort-size check: with `<10` bulk samples, neural-network methods (TAPE, Scaden) may overfit; prefer BayesPrism for small cohorts.
- Cell-type representation: cell types with `<50` reference cells produce noisy fraction estimates — collapse them into a sibling celltype before the run, or set `cellstate_key` so BayesPrism can borrow strength.
- Don't compare fractions across runs with different `max_single_cells` — the reference subsample changes and so does the signature matrix.

## Resource Map

- See [`reference.md`](reference.md) for compact copy-paste snippets per method.
- See [`references/source-grounding.md`](references/source-grounding.md) for verified `ov.bulk.Deconvolution` constructor + `.deconvolution` signatures and the four-method dispatch.
- For the opposite-direction Bulk2Single workflow (synthetic single-cell generation from bulk), see `bulk-to-single-deconvolution`.
- For spatial transcriptomics deconvolution (Tangram / cell2location / starfysh / FlashDeconv), see `spatial-tutorials`.
- For TCGA bulk preprocessing as an upstream input, see `tcga-preprocessing`.

## Examples
- "Run BayesPrism on `decov_bulk_covid_bulk()` with the bundled scRNA-seq reference; report 8-core fast-mode fractions and a stacked-bar plot."
- "Re-run with Scaden on Apple Silicon (`gpu='mps'`) and report the per-sample correlation against BayesPrism."
- "Group-average fractions by `disease` (healthy vs COVID) and plot with `ov.pl.plot_grouped_fractions`."
- "Decide which deconvolution method to pick for a 30-bulk-sample CPU-only cohort with 12 cell types in the reference."

## References
- Tutorial notebook: [`t_decov_bulk`](https://omicverse.readthedocs.io/en/latest/Tutorials-bulk/t_decov_bulk/) — COVID-19 PBMC deconvolution with BayesPrism + Scaden.
- BayesPrism: Chu *et al.* 2022, *Nature Cancer* — "Cell type and gene expression deconvolution with BayesPrism".
- Scaden: Menden *et al.* 2020, *Science Advances* — "Deep learning–based cell composition analysis from tissue expression profiles".
- TAPE: Cao *et al.* 2022, *Nature Communications* — "TAPE: an interpretable autoencoder for cell-type deconvolution".
- Live API verified — see [`references/source-grounding.md`](references/source-grounding.md).
