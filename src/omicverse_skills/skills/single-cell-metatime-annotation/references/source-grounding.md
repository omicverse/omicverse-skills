# Source Grounding — MetaTiME annotation

## Interfaces Checked

`omicverse.single.MetaTiME` and `omicverse.utils.mde`. Verified via `inspect.signature` + `inspect.getdoc` + direct reading of `omicverse/single/_metatime.py`. Cross-checked against `t_metatime.ipynb`.

## Live signatures

```python
ov.single.MetaTiME(
    adata: anndata.AnnData,
    mode: str = 'table',
)

# Methods (must run in order):
.overcluster(
    self, resolution: float = 8,
    random_state: int = 0,
    clustercol: str = 'overcluster',
)
.predictTiME(self, save_obs_name: str = 'MetaTiME')
.plot(
    self, basis: str = 'X_umap',
    cluster_key: str = 'MetaTiME',
    fontsize: int = 8,
    min_cell: int = 5,
    title=None,
    figsize: tuple = (6, 6),
    dpi: int = 80,
    frameon: bool = False,
    legend_loc=None, palette=None,
) -> (fig, ax)

ov.utils.mde(matrix: np.ndarray) -> np.ndarray
```

## Source-grounded behavior

**`MetaTiME` constructor:** stores `adata` reference (not a copy); `mode='table'` selects the bundled MeC-to-cell-state lookup table from `omicverse/external/metatime/` (the canonical mapping derived from Yi *et al.* 2023's tumor scRNA-seq atlas). `mode='cluster'` uses hard cluster-level mapping instead — rarely the right choice, kept for backward compatibility.

**`overcluster(resolution=8)`:** runs `sc.tl.leiden` at the high resolution and writes the cluster column. Resolution=8 is much higher than typical analysis resolutions because MetaTiME needs cluster-level granularity for MeC scoring; the over-clustering is *intentional*, not a regularisation choice.

**`predictTiME(save_obs_name='MetaTiME')`:**
- Computes per-cluster MeC enrichment scores against the bundled MeC matrix (~190 MeCs).
- Per cluster, picks the dominant MeC (highest score) and maps via the lookup table to a cell-state label.
- Writes:
  - `adata.obs[save_obs_name]` — fine cell-state label (e.g. 'CD8_TIL_1', 'TAM_M2_C1').
  - `adata.obs[save_obs_name + '_Major']` (canonically `'Major_MetaTiME'`) — coarse roll-up category (T_cell / Myeloid / B_cell / Tumor / Stromal / Endothelial).
- Clusters with fewer than `min_cell` cells (or where MeC scoring fails) get `'Unassigned'` in both columns.

**`plot(...)`:** custom matplotlib helper that places one label per cluster centroid using collision-aware positioning (avoids label overlap on dense embeddings). `min_cell` here filters the labels rendered, not the underlying scatter — small clusters still appear as points but without text labels.

**`ov.utils.mde`:** Minimum Distortion Embedding via `pymde` (lazily imported); produces a 2-D embedding from any input matrix. Faster than UMAP on large cohorts, qualitatively similar.

## Notebook ↔ skill alignment

| Notebook section | Skill section |
|---|---|
| Load `TiME_adata_scvi.h5ad` (scVI-corrected) | Quick Workflow §1; Input Contract |
| `ov.utils.mde(adata.obsm['X_scVI'])` → `obsm['X_mde']` | Quick Workflow §2 |
| `MetaTiME(adata, mode='table')` constructor | Quick Workflow §3 |
| `.overcluster(resolution=8, clustercol='overcluster')` | Quick Workflow §4 |
| `.predictTiME(save_obs_name='MetaTiME')` | Quick Workflow §5 |
| `.plot(cluster_key='MetaTiME', basis='X_mde', dpi=80)` | Quick Workflow §6 |
| `sc.pl.embedding(color=['Major_MetaTiME'])` | Quick Workflow §6 (alt) |

## Docstring supplementation log

`MetaTiME` (17L), `.overcluster` (10L), `.predictTiME` (11L), `.plot` (29L), `ov.utils.mde` (verified) — all already documented Numpy-style; no supplementation needed.

## Reviewer-Run Empirical Checks

- All cited functions importable: `from omicverse.single import MetaTiME; from omicverse.utils import mde` ✓
- `Major_MetaTiME` column name verified: the source writes `f'{save_obs_name}_Major'` after `predictTiME`. With default `save_obs_name='MetaTiME'`, that becomes `'Major_MetaTiME'` (matches notebook usage exactly).
- Resolution=8 is the documented and tutorial-canonical choice; other values produce poorly-segmented MeC scoring.
- The bundled MeC matrix is human-only; mouse cohorts work nominally (gene symbols overlap for many TME genes) but biological interpretation may diverge — flagged in Branch Selection.
- No live smoke run executed; the bundled `TiME_adata_scvi.h5ad` (scVI-corrected) is the canonical smoke target.
