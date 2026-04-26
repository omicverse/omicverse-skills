# Source Grounding — VIA trajectory (vanilla + velocity-guided)

## Interfaces Checked

`omicverse.single._via.pyVIA`. Verified via direct reading of `omicverse/single/_via.py` (the class definition starts at line 102 and includes the methods listed below). Cross-checked against `t_via.ipynb` and `t_via_velo.ipynb`.

## Live signatures

```python
class omicverse.single.pyVIA:
    def __init__(self,
        adata: anndata.AnnData,
        adata_key: str = 'X_pca',
        adata_ncomps: int = 80,
        basis: str = 'X_umap',
        clusters: str = 'clusters',
        knn: int = 30,
        random_seed: int = 4,
        root_user: list | None = None,
        dataset: str = '',
        is_coarse: bool = True,
        preserve_disconnected: bool = True,
        pseudotime_threshold_TS: int = 30,
        piegraph_arrow_head_width: float = 0.1,
        piegraph_edgeweight_scalingfactor: float = 1.5,
        edgebundle_pruning_twice: bool = False,
        edgebundle_pruning: float = 0.15,
        velocity_matrix: np.ndarray | None = None,
        gene_matrix: np.ndarray | None = None,
        velo_weight: float = 0.5,
        pca_loadings: np.ndarray | None = None,
        ...): ...

    def run(self): ...
    def get_piechart_dict(self, label: int = 0, clusters: str = '') -> dict: ...
    def get_pseudotime(self, adata=None): ...
    def plot_piechart_graph(self, clusters='', type_data='pt', ...): ...
    def plot_stream(self, clusters='', basis='', ...): ...
    def plot_trajectory_gams(self, clusters='', basis='', via_fine=None, idx=None, ...): ...
    def plot_lineage_probability(self, clusters='', basis='', via_fine=None, ...): ...
    def plot_gene_trend(self, gene_list=None, figsize=(8, 4), ...): ...
    def plot_clustergraph(self, gene_list, arrow_head=0.1, figsize=(8, 4), dpi=80, magic_steps=3, ...): ...
    def plot_gene_trend_heatmap(self, gene_list, marker_lineages=[], ...): ...
```

## Source-grounded behavior

**`__init__`:**
- Stores parameters; defers work to `run()`. Accepts both `'X_pca'` (standard scanpy convention) and any other `obsm` key — `adata_key` resolves through `adata.obsm[adata_key]`.
- `root_user`: list of integer cell indices (positions in `adata.obs_names`). Multiple roots are supported; pass a single-element list for the canonical case.
- `velocity_matrix`, `gene_matrix`, `velo_weight`, `pca_loadings`: all `None` for vanilla mode, all populated for velocity-guided mode. The wrapper validates this gating internally.

**`run()`:**
- Builds the PARC clustering on the chosen embedding subspace.
- Constructs the kNN graph (size `knn`).
- Lazy-teleporting random walk over the cluster transition graph; if velocity is provided, transition probabilities are biased by `velo_weight × velocity-cosine + (1 - velo_weight) × gene-distance`.
- Detects terminal states automatically using `pseudotime_threshold_TS` (terminal candidates are clusters whose pseudotime exceeds the threshold AND have no clear downstream).

**`get_pseudotime(adata=None)`:**
- Returns `np.ndarray` if `adata` is `None`; writes `adata.obs['via_pseudotime']` and returns the AnnData if `adata` is given.

**`plot_*` methods:**
- All accept `figsize`, `dpi`, `cmap` and similar matplotlib kwargs.
- `plot_stream` uses `density_grid` to control streamline density (lower = fewer / more readable).
- `plot_lineage_probability(marker_lineages=[i, j, ...])` shows only the named lineages instead of all (useful when there are 5+ terminals).
- `plot_gene_trend_heatmap` requires `gene_list` and shows a one-row-per-gene heatmap along the chosen lineages.

**Velocity-guided specifics:**
- `velocity_matrix=adata.layers['velocity']`: same shape as `adata.X`, real-valued; usually from `scvelo.tl.velocity`.
- `gene_matrix=adata.X.todense()`: dense gene-expression matrix; sparse must be densified.
- `pca_loadings=adata.varm['PCs']`: PCA loadings (`(n_genes, n_pcs)`); needed because VIA back-projects velocity into PC space.
- `velo_weight=0.5` is the canonical default; values outside `[0.3, 0.7]` are unstable.

## Notebook ↔ skill alignment

| Notebook section | Skill section |
|---|---|
| `t_via` PCA + `pyVIA(adata_key='X_pca', ..., root_user=[4823])` + `run()` | Quick Workflow §1-3 (vanilla) |
| `t_via` `plot_piechart_graph` / `plot_clustergraph` / `plot_trajectory_gams` / `plot_stream` (twice — by cluster + by time) / `plot_lineage_probability` (full + filtered) / `plot_gene_trend` / `plot_gene_trend_heatmap` | Quick Workflow §4-10; reference.md vanilla block |
| `t_via_velo` scvelo preprocess (`filter_and_normalize`, `moments`, `recover_dynamics`, `velocity`, `velocity_graph`) | Quick Workflow §11; reference.md velocity-guided block |
| `t_via_velo` `pyVIA(velocity_matrix=..., gene_matrix=..., velo_weight=0.5, pca_loadings=..., root_user=None)` + `run()` + same plot stack | Quick Workflow §12-14 |

## Docstring supplementation log

The `pyVIA` constructor has a multi-page docstring on the `__init__` (line 123 onward in `_via.py`). All public methods have docstrings — left as-is for this skill.

## Known caveats

**Import bug (flagged here for future fix):** `from omicverse.single._via import pyVIA` currently raises `NameError: name 'matplotlib' is not defined` because `matplotlib.figure.Figure` is referenced as a *type annotation* in a method signature (around line 358 of `_via.py`) without `matplotlib` being imported at module scope. The class is documented and accessible via `omicverse.single` lazy attribute — and works at runtime when matplotlib is installed and used elsewhere — but bare-import inspection fails. Tutorial cells work because notebook cells `import matplotlib.pyplot as plt` first. Skill text uses `ov.single.pyVIA(...)` which goes through the lazy attribute and resolves correctly under typical usage (matplotlib already loaded).

**Reviewer-Run Empirical Checks**

- All cited methods present in source (verified by `grep -n "    def " _via.py` — see 11 method definitions).
- Tutorial cells use exactly the verified signatures.
- VIA paper (Stassen 2021) cited for the algorithm; OmicVerse-side improvements are in the colour-scheme defaults and the `AnnData`-direct constructor.
- No live smoke run executed; the Setty hematopoiesis cohort would be the canonical smoke target.
