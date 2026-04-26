# Source Grounding — Bulk cell-type deconvolution

## Interfaces Checked

`omicverse.bulk.Deconvolution` (class + `.deconvolution` method) and the supporting demo / plotting helpers (`ov.datasets.decov_bulk_covid_*`, `ov.pl.cellproportion`, `ov.pl.plot_grouped_fractions`). Verified via `inspect.signature` + `inspect.getdoc` and direct reading of `omicverse/bulk/_deconvolution.py`. Cross-checked against `t_decov_bulk.ipynb`.

## Live signatures

```python
ov.bulk.Deconvolution(
    adata_bulk: AnnData,
    adata_single: AnnData,
    max_single_cells: int = 5000,
    celltype_key: str = 'celltype',
    cellstate_key: str | None = None,
    gpu: int | str = 0,
)

deconv_obj.deconvolution(
    method: str = 'tape',
    sep: str = '\t',
    scaler: str = 'mms',
    datatype: str = 'counts',
    genelenfile: str | None = None,
    mode: str = 'overall',
    adaptive: bool = True,
    variance_threshold: float = 0.98,
    save_model_name: str | None = None,
    batch_size: int = 128,
    epochs: int = 128,
    seed: int = 1,
    scale_size: int = 2,
    scale: bool = True,
    n_cores: int = 4,
    fast_mode: bool = True,
    pseudobulk_size: int = 2000,
    **kwargs,
) -> pd.DataFrame
```

Demo loaders + plotting helpers:
```python
ov.datasets.decov_bulk_covid_bulk()    -> AnnData    # Decov 2020 PBMC bulk
ov.datasets.decov_bulk_covid_single()  -> AnnData    # paired sc reference

ov.pl.cellproportion(adata, celltype_clusters, groupby, legend=True, ax=None)

ov.pl.plot_grouped_fractions(
    res, obs, group_key,
    color_dict=None,
    agg='mean',
    normalize=True,
    figsize=(4, 4),
)
```

## Source-grounded behavior

**Constructor:**
- Subsamples `adata_single` to `max_single_cells` if larger; the subsample is reproducible only if you set RNG state externally.
- `celltype_key` is required (the coarse annotation; output column granularity).
- `cellstate_key` is optional but used by **BayesPrism** to build a hierarchical signature (within-celltype substates pool information). Other backends ignore it.
- `gpu` semantics:
  - `int` ≥ 0: CUDA device index for Scaden / TAPE.
  - `'mps'`: Apple Silicon GPU (Scaden in the tutorial).
  - `-1`: force CPU.
  - BayesPrism is CPU-only regardless of this kwarg; tune `n_cores` instead.

**`.deconvolution(method=...)` dispatch:**
- `'tape'` (default) — TAPE autoencoder (`ov.external.tape`); accepts `mms` scaler, `adaptive=True`, `variance_threshold=0.98`. CPU-friendly.
- `'scaden'` — Scaden DNN (`ov.external.scaden`); trains on synthetic pseudo-bulk mixtures (`pseudobulk_size`); benefits from GPU; tutorial uses `scaler='ss'`, `datatype='counts'`.
- `'bayesprism'` — full Bayesian model (`ov.external.bayesprism`); `n_cores` parallel chains; `fast_mode=True` enables approximate updates. Hierarchical when `cellstate_key` is set.
- `'omicstweezer'` — joint reference-correction deconvolution (`ov.external.omicstweezer`).
- Returns a `pd.DataFrame` with samples in rows, cell types in columns; output column order is method-internal — re-order canonically with `res = res[single_ad_ref.obs[celltype_key].cat.categories]`.

**Approximate compositional constraint:** rows sum to ~1 (each sample's cell-type proportions). The methods don't strictly enforce sum-to-one; small deviations are expected.

**Built-in datasets** (`ov.datasets.decov_bulk_covid_*`):
- `decov_bulk_covid_bulk()` returns the Decov *et al.* 2020 PBMC bulk RNA-seq cohort; samples are typed by `disease` (healthy / COVID-19).
- `decov_bulk_covid_single()` is the paired single-cell reference with `cell.type.coarse` (~10 immune subsets), `cell.type.fine` (~30 sub-states), and a matching `cell.type.coarse_colors` palette in `uns`.

**`ov.pl.plot_grouped_fractions(res, obs, group_key, ...)`:** averages each cell-type column within phenotype groups (`obs[group_key]`); `agg='mean'` is default; `normalize=True` re-normalises so per-group fractions sum to 1 after averaging. Returns a stacked-bar figure.

## Notebook ↔ skill alignment

| Notebook section | Skill section |
|---|---|
| Load `decov_bulk_covid_bulk()` + `decov_bulk_covid_single()` | Quick Workflow §1 |
| UMAP coloured by `cell.type.coarse` / `cell.type.fine` / `disease` + `ov.pl.cellproportion` reference QC | Quick Workflow §2 |
| `ov.bulk.Deconvolution(... celltype_key='cell.type.coarse', cellstate_key='cell.type.fine')` | Quick Workflow §3 |
| `.deconvolution(method='bayesprism', n_cores=8, fast_mode=True)` | Quick Workflow §4; Branch Selection (BayesPrism row) |
| Re-order `res = res[single_ad_ref.obs[celltype_key].cat.categories]` | Quick Workflow §6 |
| Stacked bar with palette built from `uns['cell.type.coarse_colors']` | Quick Workflow §7 |
| Re-instantiate with `gpu='mps'` + `.deconvolution(method='scaden', scaler='ss', scale=True, datatype='counts', pseudobulk_size=2000)` | Quick Workflow §5 |

## Docstring supplementation log

- `ov.bulk.Deconvolution` (class) — 29L Numpy-style; complete.
- `.deconvolution(method=...)` — 51L; covers all four methods + every kwarg + return shape + a usage example. Complete.
- `ov.pl.plot_grouped_fractions` — 23L; covers parameter semantics. Complete.
- `ov.datasets.decov_bulk_covid_bulk` / `_single` — 3L each (one-liners). Adequate; the data origin is captured in the skill body.

No supplementation done in this skill's pass.

## Reviewer-Run Empirical Checks

- All cited functions importable: `from omicverse.bulk import Deconvolution; from omicverse.datasets import decov_bulk_covid_bulk, decov_bulk_covid_single; from omicverse.pl import cellproportion, plot_grouped_fractions` ✓
- Four `method=` values verified against the source dispatch in `omicverse/bulk/_deconvolution.py`.
- Tutorial cells use exactly the verified signatures (BayesPrism: `n_cores=8, fast_mode=True`; Scaden: `scaler='ss', scale=True, datatype='counts', pseudobulk_size=2000`).
- The COVID PBMC demo cohort with the bundled reference is the canonical smoke target; no live smoke run executed (BayesPrism's CPU run is multi-minute and Scaden's GPU run requires CUDA / MPS).
