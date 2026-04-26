---
name: omicverse-single-cell-metatime-annotation
description: Tumor microenvironment (TME) cell-state annotation via the pretrained MetaTiME meta-components (Yi et al. 2023). Three-step workflow on a batch-corrected AnnData - over-cluster at high Leiden resolution, score against pretrained MeCs, write MetaTiME / Major_MetaTiME labels into obs. Use when annotating tumor scRNA-seq cohorts, when running MetaTiME on top of an scVI / Harmony-integrated embedding, or when reproducing `t_metatime`.
---

# OmicVerse Single-Cell — MetaTiME Tumor Microenvironment Annotation

## Goal

Take a batch-corrected, dimension-reduced single-cell `AnnData` (typically from a tumor scRNA-seq cohort) and annotate it with **MetaTiME meta-components (MeCs)** — pretrained gene programs derived from millions of single cells across hundreds of tumor scRNA-seq studies. Output is a per-cell `obs['MetaTiME']` (fine-grained cell state) and `obs['Major_MetaTiME']` (coarse roll-up category). The pipeline is over-cluster → score → write.

MetaTiME's value vs. de-novo annotation is that the meta-components are **transferable across cohorts**: a TIL-1 (tissue-resident memory CD8+ T cell) annotation in your cohort means the same biology as a TIL-1 annotation in any other cohort that ran MetaTiME.

## Quick Workflow

1. Load a batch-corrected AnnData. The canonical MetaTiME demo uses scVI-corrected embedding (`adata.obsm['X_scVI']`); Harmony / Conos / scanorama corrected embeddings work equivalently.
2. (Optional) Build an MDE projection for downstream visualisation: `adata.obsm['X_mde'] = ov.utils.mde(adata.obsm['X_scVI'])`. MDE is a faster UMAP alternative bundled with OmicVerse.
3. **Construct**: `TiME_object = ov.single.MetaTiME(adata, mode='table')`. The `mode='table'` setting maps meta-component scores to cell-state labels via the bundled lookup table (default; the alternative is hard-cluster mapping but `'table'` is canonical).
4. **Over-cluster at high resolution**: `TiME_object.overcluster(resolution=8, clustercol='overcluster')`. Writes `adata.obs['overcluster']`. Resolution=8 (much higher than typical analysis resolutions of 0.5–1.5) is required — MetaTiME's MeC scoring is computed at the cluster level, so over-clustering preserves cell-state granularity.
5. **Predict**: `TiME_object.predictTiME(save_obs_name='MetaTiME')`. Scores each cluster against all MeCs, assigns the dominant MeC, and writes both fine-grained (`obs['MetaTiME']`) and major-category (`obs['Major_MetaTiME']`) labels.
6. **Visualise**: `TiME_object.plot(cluster_key='MetaTiME', basis='X_mde', dpi=80)` — built-in plot with collision-aware label placement, OR fall back to `sc.pl.embedding(adata, basis='X_mde', color=['Major_MetaTiME'], frameon=False)` for a standard scanpy figure.

## Interface Summary

```python
ov.single.MetaTiME(
    adata: AnnData,
    mode: str = 'table',           # 'table' (canonical lookup) or 'cluster' (hard cluster mapping)
)
# Constructor; stores adata reference. No work done yet.
```

Methods (must run in order):
- `TiME_object.overcluster(resolution: float = 8, random_state: int = 0, clustercol: str = 'overcluster')` — high-resolution Leiden clustering. Writes `adata.obs[clustercol]`.
- `TiME_object.predictTiME(save_obs_name: str = 'MetaTiME')` — scores clusters against MeCs, writes `adata.obs[save_obs_name]` (fine) and `adata.obs[save_obs_name + '_Major']` (typically `'Major_MetaTiME'`).
- `TiME_object.plot(basis: str = 'X_umap', cluster_key: str = 'MetaTiME', fontsize: int = 8, min_cell: int = 5, title=None, figsize: tuple = (6, 6), dpi: int = 80, frameon: bool = False, legend_loc=None, palette=None) → (fig, ax)` — collision-adjusted labelled scatter on the chosen embedding.

Helper:
- `ov.utils.mde(matrix) → ndarray` — Minimum Distortion Embedding; faster than UMAP, similar quality. Pass any low-D representation (PCA / scVI / Harmony / scanorama) to project to 2D.

## Boundary

**Inside scope:**
- Three-step MetaTiME annotation on a batch-corrected `AnnData`.
- Over-clustering + MeC scoring + label writeback.
- Built-in plot with collision-aware labels.
- Major-category roll-up (Tumor / T cell / Myeloid / Stromal / etc.).

**Outside scope — separate skill:**
- Batch correction itself (scVI / Harmony / Conos / scanorama) — see existing batch-integration skill.
- Other annotation tools (CellTypist, gpt4celltype, scsa) — see existing `single-cell-annotation` skill.
- Cell Ontology mapping of free-text annotations — see `omicverse-single-cell-cellmatch-ontology`.
- Multi-annotator consensus voting — see `omicverse-single-cell-cellvote-consensus`.
- Non-tumor cohorts: MetaTiME's MeCs are derived from tumor data and are most useful in TME context. On normal-tissue cohorts, MeC labels still resolve but the biological interpretation is less direct.

## Branch Selection

**`mode='table'` vs `mode='cluster'`**
- `'table'` (default) — uses the bundled lookup table to map MeC scores to cell-state labels. Canonical and stable across MetaTiME versions.
- `'cluster'` — hard cluster-level mapping based on majority MeC; less granular, rarely needed.

**`resolution` for `overcluster`**
- 8 (tutorial default) — appropriate for cohorts of 5k–50k cells.
- 10–12 — for very large cohorts (>100k cells); preserves more sub-states.
- 4–6 — for small cohorts (<5k cells); avoids over-fragmentation.
- Below 4 — too coarse; MeC scoring loses resolution.

**`basis` for `plot`**
- `'X_mde'` (recommended) — fast and visually similar to UMAP.
- `'X_umap'` — standard.
- `'X_tsne'` — when the cohort has clear cluster separation.
- `'X_scVI'` (raw 30-D scVI) — won't work; needs a 2-D embedding.

**`fontsize` and `min_cell` for `plot`**
- `fontsize=8` (default) — readable on most figures.
- `min_cell=5` — only label clusters with ≥5 cells. Increase to 20–50 for larger cohorts to reduce label clutter.

**`Major_MetaTiME` vs `MetaTiME`**
- `Major_MetaTiME` — coarse categories (e.g. T_cell, Myeloid, Tumor, Stromal). Use for high-level cohort summaries.
- `MetaTiME` — fine cell states (e.g. CD8_TIL_1, TAM_M2_C1, Treg_TIM3+). Use for biological interpretation.
- Always report both in figures: major labels as the colour, fine labels as text.

## Input Contract

- `AnnData` with a batch-corrected low-D representation in `obsm` (`'X_scVI'`, `'X_harmony'`, `'X_pca_corrected'`, etc.). MetaTiME expects integrated data — running on a single-batch cohort works but skip step 2.
- Counts in `.X` should be log-normalised (`sc.pp.normalize_total` + `sc.pp.log1p`); MeC scoring is on log-CPM scale.
- Gene symbols (not IDs) in `var_names` — the MeC lookup is keyed by symbol.
- Cohort must be human (the bundled MeCs are human-derived); mouse cohorts need symbol-translation or a different MeC set (not bundled).

## Minimal Execution Patterns

```python
import omicverse as ov
import scanpy as sc

ov.utils.ov_plot_set()

# 1) Load batch-corrected AnnData (scVI corrected here)
adata = sc.read('TiME_adata_scvi.h5ad')

# 2) MDE projection for visualisation
adata.obsm['X_mde'] = ov.utils.mde(adata.obsm['X_scVI'])

# 3) Construct + over-cluster + predict
TiME_object = ov.single.MetaTiME(adata, mode='table')
TiME_object.overcluster(resolution=8, clustercol='overcluster')
TiME_object.predictTiME(save_obs_name='MetaTiME')

# obs now has:  'overcluster', 'MetaTiME', 'Major_MetaTiME'
print(adata.obs[['MetaTiME', 'Major_MetaTiME']].value_counts().head(10))

# 4) Visualise — built-in collision-adjusted labels
fig, ax = TiME_object.plot(
    cluster_key='MetaTiME',
    basis='X_mde',
    dpi=80,
)

# Or a standard scanpy embedding plot for the major categories
sc.pl.embedding(
    adata, basis='X_mde',
    color=['Major_MetaTiME'],
    frameon=False, ncols=1,
)
```

## Validation

- After `overcluster(resolution=8)`: `adata.obs['overcluster'].nunique()` should be in the dozens–low hundreds. If it's <10, the resolution is too low for MetaTiME (over-clustering is the *point*).
- After `predictTiME`: both `MetaTiME` and `Major_MetaTiME` columns are populated, no nulls. If many cells get `'Unassigned'`, the MeC scoring failed — check that gene symbols match human convention and that `.X` is log-normalised.
- Sanity check: `Major_MetaTiME` should include the canonical TME categories you expect (Tumor, T_cell, Myeloid, B/Plasma, Stromal, Endothelial, etc.). If the cohort is tumor-only and only `Tumor` appears as a major, the MeC scoring is wrong (probably input is scaled instead of log).
- The `MetaTiME` plot is busy by default — clusters smaller than `min_cell` are unlabelled to reduce clutter.
- Don't run MetaTiME on a non-batch-corrected cohort — batch effects dominate the over-clustering and produce a per-batch fragmentation rather than per-cell-state.

## Resource Map

- See [`reference.md`](reference.md) for compact copy-paste snippets.
- See [`references/source-grounding.md`](references/source-grounding.md) for the verified `MetaTiME` constructor + method signatures and the bundled MeC lookup behavior.
- For batch correction (scVI / Harmony / Conos / scanorama) before MetaTiME, see the existing `single-cell-batch-integration` skill.
- For non-MetaTiME annotators (CellTypist / gpt4celltype / scsa), see the existing `single-cell-annotation` skill.
- For Cell Ontology mapping of MetaTiME labels, see `omicverse-single-cell-cellmatch-ontology`.

## Examples
- "Annotate a scVI-corrected tumor scRNA-seq AnnData with MetaTiME at over-cluster resolution 8 and plot fine and major categories on MDE."
- "Compare MetaTiME's `Major_MetaTiME` rollup with manual marker-based annotation; flag clusters where they disagree."
- "Run MetaTiME on a 100k-cell pan-cancer cohort with `resolution=10` and report the cell-state distribution."
- "Diagnose why MetaTiME assigned all cells to a single major category — likely the input wasn't log-normalised."

## References
- Tutorial notebook: [`t_metatime.ipynb`](https://omicverse.readthedocs.io/en/latest/Tutorials-single/t_metatime/) — TIL annotation walkthrough on a scVI-corrected cohort.
- MetaTiME paper: Yi *et al.* 2023 — "MetaTiME integrates single-cell gene expression to characterize the meta-components of the tumor immune microenvironment".
- Live API verified — see [`references/source-grounding.md`](references/source-grounding.md).
