---
name: omicverse-single-cell-via-trajectory
description: VIA (Stassen 2021) single-cell trajectory inference via `omicverse.single.pyVIA` - PARC-clustering + lazy-teleporting random-walk + automated terminal-state detection + temporal gene-trend GAMs. Two modes - vanilla (gene-distance only) and RNA-velocity-guided (`velocity_matrix`, `gene_matrix`, `velo_weight`). Use when running VIA / scVelo+VIA on AnnData, when reproducing `t_via` or `t_via_velo`, or when picking the right `velo_weight` for velocity-driven topology.
---

# OmicVerse Single-Cell — VIA Trajectory (with optional RNA velocity)

## Goal

Take a preprocessed annotated single-cell `AnnData` and run **VIA** (Stassen *et al.* 2021, Nature Communications) for: (a) topology construction, (b) pseudotime, (c) automated terminal-state detection, (d) temporal gene-trend visualisation along inferred lineages. VIA's signature feature vs. classical trajectory methods is **automated terminal-state prediction** — you don't have to specify endpoints.

Two modes:
1. **Vanilla VIA** (`t_via`) — gene-distance only; root cell required.
2. **Velocity-guided VIA** (`t_via_velo`) — pass `velocity_matrix=adata.layers['velocity']`, `gene_matrix=adata.X.todense()`, `velo_weight=0.5`. Velocity guides the random-walk transition probabilities; root cell can be left automatic (`root_user=None`).

`pyVIA` here is OmicVerse's wrapper — it improves the original VIA's colour scheme to default OmicVerse-friendly choices and lets the user pass `AnnData` directly instead of explicit numpy arrays.

## Quick Workflow

### Vanilla VIA (`t_via`)

1. Load preprocessed AnnData with cell-type labels and a low-dim representation (`X_pca` recommended; pre-compute with `sc.tl.pca(adata, svd_solver='arpack', n_comps=200)` if not present).
2. Construct: `v0 = ov.single.pyVIA(adata=adata, adata_key='X_pca', adata_ncomps=80, basis='tsne', clusters='label', knn=30, random_seed=4, root_user=[<cell_index>])`.
3. Run: `v0.run()`. Internally builds the PARC clustering, the lazy-teleporting random-walk, and the lineage probability matrix.
4. **Topology view**: `v0.plot_piechart_graph(clusters='label', cmap='Reds', show_legend=False, ax_text=False, fontsize=4)` — pie-chart graph of cluster transitions.
5. **GAM-fitted trajectory curves**: `v0.plot_trajectory_gams(basis='tsne', clusters='label', draw_all_curves=False)` — major lineage curves overlaid on the embedding.
6. **Streamplot**: `v0.plot_stream(basis='tsne', clusters='label', density_grid=0.8, scatter_size=30, ...)` — directional flow on the embedding.
7. **Lineage probabilities**: `v0.plot_lineage_probability(figsize=(8, 4))` for all detected lineages, or `marker_lineages=[2, 3]` for a subset.
8. **Per-gene trends along lineages**: `v0.plot_gene_trend(gene_list=marker_genes, figsize=(8, 6))` — GAM curves per gene, one panel per lineage.
9. **Per-lineage gene heatmap**: `v0.plot_gene_trend_heatmap(gene_list=marker_genes, marker_lineages=[2], figsize=(4, 4))`.
10. **Cluster-graph with gene overlays**: `v0.plot_clustergraph(gene_list=top_genes[:4], figsize=(12, 3))` — gene expression coloured on the cluster transition graph.

### Velocity-guided VIA (`t_via_velo`)

11. Pre-compute RNA velocity with `scvelo` (or `dynamo` etc.) so `adata.layers['velocity']` is populated; PCA must also be computed (and `adata.varm['PCs']` available).
12. Construct: `v0 = ov.single.pyVIA(adata, adata_key='X_pca', adata_ncomps=n_pcs, basis='X_umap', clusters='clusters', knn=20, root_user=None, is_coarse=True, preserve_disconnected=True, pseudotime_threshold_TS=50, piegraph_arrow_head_width=0.15, piegraph_edgeweight_scalingfactor=2.5, velocity_matrix=adata.layers['velocity'], gene_matrix=adata.X.todense(), velo_weight=0.5, edgebundle_pruning_twice=False, edgebundle_pruning=0.15, pca_loadings=adata.varm['PCs'], random_seed=42)`.
13. Run: `v0.run()`. Velocity now guides the transition probabilities; `root_user=None` is OK because terminal states (and implicitly the root) are inferred.
14. Same plots as vanilla (`plot_piechart_graph`, `plot_trajectory_gams`, `plot_stream`, `plot_lineage_probability`).

## Interface Summary

```python
ov.single.pyVIA(
    adata: AnnData,
    adata_key: str = 'X_pca',
    adata_ncomps: int = 80,
    basis: str = 'X_umap',
    clusters: str = 'clusters',
    knn: int = 30,
    random_seed: int = 4,
    root_user: list | None = None,           # cell indices to seed VIA's root
    dataset: str = '',
    is_coarse: bool = True,
    preserve_disconnected: bool = True,
    pseudotime_threshold_TS: int = 30,       # bumped to 50 in t_via_velo
    piegraph_arrow_head_width: float = 0.1,
    piegraph_edgeweight_scalingfactor: float = 1.5,
    edgebundle_pruning_twice: bool = False,
    edgebundle_pruning: float = 0.15,
    # velocity-guided mode (all None for vanilla VIA)
    velocity_matrix: np.ndarray | None = None,    # adata.layers['velocity']
    gene_matrix: np.ndarray | None = None,        # adata.X.todense()
    velo_weight: float = 0.5,                     # 0 = ignore velocity, 1 = pure velocity
    pca_loadings: np.ndarray | None = None,       # adata.varm['PCs']
    ...
)
```

Methods (after `v0.run()`):
- `v0.run()` — fit; populates internal state.
- `v0.get_pseudotime(adata=None) → np.ndarray | adata` — per-cell pseudotime; can write into `adata.obs['via_pseudotime']` when given.
- `v0.get_piechart_dict(label=0, clusters='') → dict` — per-cluster lineage probabilities, the data backing `plot_piechart_graph`.
- `v0.plot_piechart_graph(clusters='', type_data='pt', ...) → (fig, ax, ax1)` — pie-chart cluster transition graph.
- `v0.plot_clustergraph(gene_list, arrow_head=0.1, figsize=(8, 4), dpi=80, magic_steps=3, ...)` — cluster graph with per-gene expression overlay.
- `v0.plot_stream(clusters='', basis='', density_grid=0.8, scatter_size=30, scatter_alpha=0.3, linewidth=0.5, color_scheme='time', min_mass=1, cutoff_perc=5, marker_edgewidth=0.1, density_stream=2, smooth_transition=1, smooth_grid=0.5, ...)` — directional flow streamplot.
- `v0.plot_trajectory_gams(clusters='', basis='', via_fine=None, idx=None, draw_all_curves=False, ...)` — GAM-fitted lineage curves.
- `v0.plot_lineage_probability(clusters='', basis='', via_fine=None, marker_lineages=None, figsize=(8, 4), ...)`.
- `v0.plot_gene_trend(gene_list, figsize=(8, 4), ...)` — GAM fits per gene per lineage.
- `v0.plot_gene_trend_heatmap(gene_list, marker_lineages=[], ...)`.

Dataset helpers:
- `ov.single.scRNA_hematopoiesis()` — Setty 2019 hematopoiesis demo (one of two canonical VIA test cohorts).

## Boundary

**Inside scope:**
- VIA fitting on `AnnData` with `pyVIA` constructor + `.run()`.
- Vanilla and velocity-guided modes.
- All VIA plotting helpers (piechart, stream, trajectory_gams, lineage_probability, gene_trend, gene_trend_heatmap, clustergraph).

**Outside scope — separate skill:**
- Other trajectory inference methods (Monocle2 / diffusion / slingshot / palantir / scTour) — see those skills.
- RNA velocity computation itself (scvelo / dynamo / latentvelo / graphvelo) — see `single-cell-rna-velocity` (existing skill).
- CellRank fate maps from RNA velocity (a different post-velocity downstream) — see `omicverse-single-cell-cellrank-fate`.
- Multi-omic VIA (StaVIA on metabolomics / RNA — that's covered by `single-multiomics` skill).

## Branch Selection

**Vanilla vs velocity-guided**
- Vanilla (`velocity_matrix=None`, `root_user=[<cell>]`): fastest; works on any AnnData with a sane PCA. The `root_user` cell must be picked from the earliest progenitor — otherwise pseudotime is reversed.
- Velocity-guided (`velocity_matrix=...`, `root_user=None`): more accurate when velocity is reliable; **automatic root inference**. Use when the cohort has good splice-junction coverage (Smart-seq2 or similar) and `scvelo.tl.recover_dynamics` ran successfully. Don't use on shallow 10x runs with sparse spliced/unspliced separation.

**`velo_weight` (velocity-guided mode only)**
- 0.0 — equivalent to vanilla (velocity ignored).
- 0.5 — tutorial default; balances velocity and gene-distance information.
- 1.0 — pure velocity. Often *too* noisy on shallow runs because of poorly-estimated velocities for low-expression genes.
- The Stassen paper recommends 0.3–0.7 as the safe range. Below 0.3 the velocity is barely contributing; above 0.7 noisy velocity estimates can drive bad transitions.

**`adata_ncomps` (PCs to use)**
- 80 (vanilla default for ~5k-cell hematopoietic): captures dominant lineage structure.
- For larger / more heterogeneous cohorts: 100–200.
- For very small cohorts (<1k cells): 30–50.
- Stick close to 50–100 for most studies; pushing higher rarely helps and slows the random walk.

**`knn`**
- 30 — vanilla default; appropriate for cohorts ~1k–10k cells.
- 20 — velocity-guided tutorial default (the velocity is already smoothing transitions, so a tighter neighbour graph is fine).
- 50+ — for cohorts >50k cells; larger graph stabilises the cluster transition probabilities at the cost of sensitivity to rare branches.

**`random_seed`** — set explicitly. VIA's PARC clustering and random walk both have stochastic components.

**`root_user`** — pass a cell *index* (integer position in `adata.obs_names`) of an early progenitor. For the Setty hematopoiesis demo, `[4823]` is the canonical CD34+ HSC. With velocity-guided mode, leave `None`.

**`is_coarse=True`** — coarse-grained (faster, ~tens of clusters); leave `True` unless you specifically need cell-resolution VIA (rare).

**`preserve_disconnected=True`** — keep disconnected cluster components (e.g., contaminating cell types). Set `False` to force a single connected component (rarely the right choice for real data).

**`pseudotime_threshold_TS`** (terminal-state pseudotime threshold)
- 30 — vanilla default; aggressive.
- 50 — velocity-guided tutorial default; less aggressive (more candidate terminals).
- Lower → more terminal states detected (some may be spurious).
- Higher → fewer, more confident terminals.

**`edgebundle_pruning`** controls how many transition edges survive in the final piechart graph. 0.15 (default) keeps a moderate density. Lower → fewer edges (cleaner figure); higher → more (more biology, busier figure).

**Plotting `color_scheme`**
- `'time'` (default for `plot_stream`) — colour by pseudotime; flow goes from cool to warm.
- Pass `clusters='label'` to colour by cell-type cluster instead.

## Input Contract

- `AnnData` with `obs[<celltype>]` populated; PCA in `obsm['X_pca']` (or whatever `adata_key` says); 2-D embedding in `obsm[basis]` (`tsne` / `umap` / etc.).
- `var_names` are gene symbols; gene-trend / gene-heatmap / clustergraph plots index by gene name.
- For velocity-guided mode: `adata.layers['velocity']` populated by scvelo (or compatible); `adata.varm['PCs']` (PCA loadings) — write with `sc.tl.pca` first.
- `root_user`: integer index, not a string cell name.
- For sparse `.X`, the velocity-guided path may need `adata.X.todense()` because `gene_matrix` is consumed as a dense numpy array.

## Minimal Execution Patterns

```python
# Vanilla VIA
import omicverse as ov
import scanpy as sc
import matplotlib.pyplot as plt

ov.utils.ov_plot_set()

adata = ov.single.scRNA_hematopoiesis()
sc.tl.pca(adata, svd_solver='arpack', n_comps=200)

v0 = ov.single.pyVIA(
    adata=adata,
    adata_key='X_pca', adata_ncomps=80,
    basis='tsne',
    clusters='label',
    knn=30, random_seed=4,
    root_user=[4823],
)
v0.run()

# Diagnostics
v0.plot_piechart_graph(clusters='label', cmap='Reds', show_legend=False,
                       ax_text=False, fontsize=4)
v0.plot_trajectory_gams(basis='tsne', clusters='label', draw_all_curves=False)
v0.plot_stream(basis='tsne', clusters='label', density_grid=0.8,
               scatter_size=30, scatter_alpha=0.3, linewidth=0.5)
v0.plot_lineage_probability(figsize=(8, 4))

# Per-gene trends
markers = ['IL3RA', 'IRF8', 'GATA1', 'GATA2', 'ITGA2B', 'MPO',
           'CD79B', 'SPI1', 'CD34', 'CSF1R', 'ITGAX']
v0.plot_gene_trend(gene_list=markers, figsize=(8, 6))
v0.plot_gene_trend_heatmap(gene_list=markers, marker_lineages=[2], figsize=(4, 4))

# Cluster graph with gene overlays
v0.plot_clustergraph(gene_list=markers[:4], figsize=(12, 3))
```

```python
# Velocity-guided VIA
import scvelo as scv
import scanpy as sc

# Pre-compute RNA velocity (separate skill)
scv.pp.filter_and_normalize(adata, min_shared_counts=10, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)
sc.tl.pca(adata, n_comps=30)
n_pcs = 30
sc.tl.umap(adata)

# Construct VIA with velocity
v0 = ov.single.pyVIA(
    adata=adata,
    adata_key='X_pca', adata_ncomps=n_pcs,
    basis='X_umap',
    clusters='clusters',
    knn=20, root_user=None,                  # automatic root
    is_coarse=True, preserve_disconnected=True,
    pseudotime_threshold_TS=50,
    piegraph_arrow_head_width=0.15,
    piegraph_edgeweight_scalingfactor=2.5,
    velocity_matrix=adata.layers['velocity'],
    gene_matrix=adata.X.todense(),
    velo_weight=0.5,
    edgebundle_pruning_twice=False,
    edgebundle_pruning=0.15,
    pca_loadings=adata.varm['PCs'],
    random_seed=42,
)
v0.run()

v0.plot_piechart_graph(clusters='clusters', cmap='Reds', show_legend=False,
                       ax_text=False, fontsize=4)
v0.plot_trajectory_gams(basis='X_umap', clusters='clusters', draw_all_curves=False)
v0.plot_stream(basis='X_umap', clusters='clusters', density_grid=0.8,
               scatter_size=30, scatter_alpha=0.3, linewidth=0.5)
v0.plot_lineage_probability()
```

## Validation

- After `v0.run()`: VIA prints a summary listing detected terminal states. The number of terminals should be biologically sensible (2–6 for hematopoiesis; >10 means `pseudotime_threshold_TS` is too low).
- `v0.get_pseudotime()` returns finite values for every cell. NaN values indicate disconnected cluster components — set `preserve_disconnected=False` (with caveats) or revisit the kNN graph.
- `plot_piechart_graph` should show a clear arrow-direction structure: arrows from progenitor cluster outward. If arrows are random, the root cell is wrong.
- For velocity-guided mode: compare the VIA topology against `scvelo.pl.velocity_embedding_stream` — the directional flows should agree on dominant lineages. Disagreement implies one of (a) velocity mis-estimated, (b) `velo_weight` wrong, (c) PCA loadings stale.
- Per-gene trends: marker genes should peak at biologically expected lineages (e.g. GATA1 → erythroid, MPO → myeloid). If markers don't peak in the expected lineage, the VIA terminals are mis-assigned and the lineage indices in `plot_gene_trend` need to be reinterpreted.
- Disagreement between vanilla and velocity-guided VIA on the same cohort: signals that velocity is adding (or hiding) signal — investigate via `scvelo.pl.velocity_confidence`.

## Resource Map

- See [`reference.md`](reference.md) for compact copy-paste snippets for both modes.
- See [`references/source-grounding.md`](references/source-grounding.md) for verified `pyVIA` constructor + method signatures and a known-import-bug caveat.
- For RNA velocity computation, see the existing `single-cell-rna-velocity` skill.
- For other trajectory backends (Monocle2 / diffusion / slingshot / palantir / scTour), see the respective per-method skills.
- For CellRank fate maps from RNA velocity, see `omicverse-single-cell-cellrank-fate`.

## Examples
- "Run vanilla VIA on `ov.single.scRNA_hematopoiesis()` with `root_user=[4823]` and plot the pie-chart graph + lineage probabilities + gene trends for GATA1 / MPO / SPI1."
- "Run velocity-guided VIA with `velo_weight=0.5` on a scvelo-processed AnnData; verify the terminal states match the dominant lineages in the velocity stream."
- "Diagnose why VIA pseudotime is reversed — likely the `root_user` cell is in a terminal cluster instead of a progenitor."
- "Compare vanilla and velocity-guided VIA on the same cohort; identify lineages where they disagree."

## References
- Tutorial notebooks:
  - [`t_via.ipynb`](https://omicverse.readthedocs.io/en/latest/Tutorials-single/t_via/) — vanilla VIA on Setty hematopoiesis.
  - [`t_via_velo.ipynb`](https://omicverse.readthedocs.io/en/latest/Tutorials-single/t_via_velo/) — velocity-guided VIA on the dentate-gyrus / pancreas style cohort.
- VIA paper: Stassen *et al.* 2021, *Nature Communications* — "Generalized and scalable trajectory inference in single-cell omics data with VIA".
- Live API verified — see [`references/source-grounding.md`](references/source-grounding.md).
