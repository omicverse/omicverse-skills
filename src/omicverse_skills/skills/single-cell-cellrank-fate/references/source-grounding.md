# Source Grounding — CellRank fate maps

## Interfaces Checked

This skill is a **post-velocity** workflow whose core algorithm lives in the external `cellrank` PyPI package (Lange *et al.* 2022). The OmicVerse-side contribution is the visualisation stack (`ov.pl.branch_streamplot`, `ov.pl.dynamic_heatmap`) and the GAM-fitting backend (`ov.single.dynamic_features`) — those were verified via `inspect.signature` + `inspect.getdoc` and direct reading of `omicverse/single/dynamic_features.py` + `omicverse/pl/_branch_streamplot.py` + `omicverse/pl/_dynamic_*.py`. CellRank API is documented from notebook usage cross-checked against CellRank 2's published API.

## Live signatures — CellRank (external, `pip install cellrank`)

```python
cr.kernels.VelocityKernel(adata).compute_transition_matrix() -> Kernel
cr.kernels.ConnectivityKernel(adata).compute_transition_matrix() -> Kernel
# Kernel arithmetic via operator overloads:
combined_kernel = 0.8 * vk + 0.2 * ck

cr.estimators.GPCCA(kernel)
    .compute_macrostates(n_states=N, cluster_key='clusters')
    .predict_terminal_states(n_states=k)
    .fate_probabilities  # property -> pd.DataFrame (n_cells, n_terminals)
    .macrostates         # property -> pd.Categorical
    .terminal_states     # property -> pd.Categorical
```

## Live signatures — OmicVerse-side (`ov.pl` + `ov.single`)

```python
ov.pl.branch_streamplot(
    adata, *,
    group_key: str,
    pseudotime_key: str,
    trunk_groups: list[str] = None,
    branch_center: float = 0.5,
    figsize: tuple = ...,
    xlabel: str = None,
    show: bool = True,
    ...,
) -> (fig, ax)

ov.single.dynamic_features(
    adata, genes, pseudotime: str, *,
    layer: str | None = None,
    groupby: str | None = None,
    groups: Sequence[str] | None = None,
    distribution: str = 'normal',
    link: str = 'identity',
    n_splines: int = 8,
    spline_order: int = 3,
    grid_size: int = 200,
    confidence_level: float = 0.95,
    min_cells: int = 20,
    min_variance: float = 1e-08,
    store_raw: bool = False,
    raw_obs_keys: list[str] | None = None,
    key_added: str | None = 'dynamic_features',
    verbose: bool = True,
) -> DynamicFeaturesResult

ov.pl.dynamic_trends(
    res, *,
    genes,
    compare_features: bool = False,
    compare_groups: bool = False,
    split_time: float | None = None,
    shared_trunk: bool = True,
    add_point: bool = True,
    point_color_by: str | None = None,
    line_style_by: str | None = None,
    figsize, linewidth, ncols, legend_loc, legend_fontsize, title, ...
)

ov.pl.dynamic_heatmap(
    adata, *,
    pseudotime: str,
    var_names: list | dict,         # dict triggers multi-module (panel) layout
    cell_annotation: str | None = None,
    layer: str | None = None,
    use_fitted: bool = True,
    cell_bins: int = 200,
    smooth_window: int = ...,
    fitted_window: int = ...,
    standard_scale: str | None = 'var',
    cmap: str = 'RdBu_r',
    order_by: str = 'peak',
    show_row_names: bool = True,
    figsize: tuple = ...,
    border: bool = False,
    show: bool = True,
)

ov.pl.embedding(adata, *, basis='X_umap', color, cmap, frameon, ...)
```

## Source-grounded behavior

**Kernel mix:** CellRank's `Kernel` class supports `__add__` / `__rmul__`, so `0.8 * vk + 0.2 * ck` returns a new combined kernel whose `transition_matrix` is the convex combination. Weights normalise to 1 internally; passing `1.5 * vk + 0.5 * ck` is silently re-normalised to `0.75 * vk + 0.25 * ck`.

**GPCCA macrostates:** spectral decomposition on the combined kernel; `n_states` is the number of macrostates to retain, NOT the number of terminals. `predict_terminal_states(n_states=k)` then selects `k` candidate terminals from the macrostates (excludes trunks based on out-flow signatures).

**`fate_probabilities`** is `(n_cells, n_terminals)`; column names are the terminal-state labels (pulled from `g.terminal_states.cat.categories`). Pulling `g.fate_probabilities['Beta']` returns a column-vector ndarray-like; `np.asarray(...).ravel()` flattens it to `(n_cells,)`.

**`branch_streamplot`** in OmicVerse: combines a streamplot (with trunk + branch directionality from `pseudotime_key`) with cluster colour overlay. `trunk_groups` lists clusters that should be drawn as a single trunk; `branch_center` is the x-position (in normalised pseudotime [0, 1]) where the trunk splits.

**`dynamic_features(layer='Ms')`:** the `layer` kwarg routes the GAM fit to `adata.layers[layer]` rather than `adata.X`. For velocity workflows this MUST be `'Ms'` (smoothed spliced from scvelo) — using `adata.X` (typically log1p-normalised) gives systematically wrong scales; using `'spliced'` gives noisy raw counts.

**`dynamic_heatmap(var_names: dict)`:** when `var_names` is a `dict[str, list[str]]`, the heatmap renders one row-group per dict entry with the entry name as a panel label. This is the multi-module pancreas figure pattern.

## Notebook ↔ skill alignment

| Notebook section | Skill section |
|---|---|
| Kernel construction (`vk + ck`) | Quick Workflow §2; Branch Selection (kernel mix) |
| `GPCCA.compute_macrostates(n_states=6)` + `predict_terminal_states(n_states=4)` | Quick Workflow §3-4 |
| Per-cell fate probability → `obs['beta_fate']` + UMAP plot | Quick Workflow §5; reference.md fate-probability block |
| `ov.pl.branch_streamplot(trunk_groups, branch_center=0.62)` | Quick Workflow §6 |
| `dynamic_features(adata_beta, layer='Ms')` + `dynamic_trends(compare_features=True)` | Quick Workflow §7-8 |
| Branch-aware Alpha-vs-Beta with `shared_trunk=True` | Quick Workflow §9; reference.md |
| Multi-module `dynamic_heatmap(var_names=modules)` | Quick Workflow §10 |

## Docstring supplementation log

`ov.single.dynamic_features` (80L) — already well-documented. `ov.pl.branch_streamplot`, `ov.pl.dynamic_trends`, `ov.pl.dynamic_heatmap` are documented in their respective `omicverse/pl/*.py` files. No supplementation needed for this skill.

## Reviewer-Run Empirical Checks

- All cited `ov` functions importable: `from omicverse.single import dynamic_features; from omicverse.pl import branch_streamplot, dynamic_trends, dynamic_heatmap, embedding` ✓
- CellRank API conforms to the published 2.0 interface; the tutorial's exact calls (`VelocityKernel`, `ConnectivityKernel`, `0.8 * vk + 0.2 * ck`, `GPCCA`, `compute_macrostates`, `predict_terminal_states`, `fate_probabilities`) are the canonical CellRank 2 examples.
- `layer='Ms'` requirement verified against scvelo convention (smoothed spliced layer is what the GAM should consume post-velocity).
- No live smoke run executed; the pancreas endocrinogenesis cohort (`scv.datasets.pancreas()`) is the canonical smoke target for both velocity computation and CellRank fate mapping.

## Caveats

- CellRank is an external dependency (`pip install cellrank`); the skill description and Input Contract make this explicit.
- The `Ms` layer naming is scvelo-specific; if velocity was computed with `dynamo`, the equivalent layer is `Ms_dyn` (or similar) — adjust `dynamic_features(layer=...)` accordingly. Not all velocity backends produce a smoothed-spliced layer; `dynamic_features` will fall back to `adata.X` if `layer` is `None`.
