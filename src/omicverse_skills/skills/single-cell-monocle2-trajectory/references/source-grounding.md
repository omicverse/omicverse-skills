# Source Grounding — Monocle2 trajectory

## Interfaces Checked

`omicverse.single.Monocle` and the trajectory-related plotting / GAM API. Verified via `inspect.signature` + `inspect.getdoc` + direct reading of `omicverse/single/_monocle*.py` (or wherever the Monocle implementation lives) and `omicverse/single/dynamic_features.py`. Cross-checked against `t_traj_monocle2.ipynb` (renamed from `t_traj_monocle2_olsson` in tutorials commit `12ce9f4`).

## Live signatures

```python
ov.single.Monocle(adata: AnnData)
# class with .adata attribute (a copy of the input)

# Methods (verified from notebook usage; full source in omicverse/single/):
mono.preprocess(...)
mono.select_ordering_genes(max_genes: int = 1000)
mono.plot_ordering_genes(figsize: tuple = ...)
mono.fit_trajectory(...)         # invoked implicitly by plot_trajectory
mono.plot_trajectory(
    color_by='State',
    cell_size: float = 2.4,
    cell_link_size: float = 0.5,
    show_branch_points: bool = False,
    figsize: tuple = (6, 4),
    palette=None,
    ...
)
mono.differential_gene_test(cores: int = -1) -> pd.DataFrame
# columns: pval, qval, status ('OK'|'FAIL'), gene_id (index)
mono.BEAM(branch_point: int = 1, cores: int = -1) -> pd.DataFrame
# same column schema as differential_gene_test
```

```python
ov.single.dynamic_features(
    data: AnnData | Mapping[str, AnnData],
    genes: Sequence[str] | str | None = None,
    pseudotime: str = 'pseudotime',
    *,
    groupby: str | Mapping[str, str] | None = None,
    groups: Sequence[str] | Mapping[str, Sequence[str]] | None = None,
    layer: str | None = None,
    use_raw: bool = False,
    subsets: Sequence[str] | Mapping[str, Sequence[str]] | None = None,
    weights: str | Sequence[float] | Mapping[str, Any] | None = None,
    distribution: str = 'normal',
    link: str = 'identity',
    n_splines: int = 8,
    spline_order: int = 3,
    grid_size: int = 200,
    confidence_level: float = 0.95,
    min_cells: int = 20,
    min_variance: float = 1e-08,
    store_raw: bool = False,
    raw_obs_keys: Sequence[str] | str | Mapping[str, Sequence[str] | str] | None = None,
    key_added: str | None = 'dynamic_features',
    verbose: bool = True,
) -> DynamicFeaturesResult

ov.pl.branch_streamplot(adata, *, group_key, pseudotime_key, show=True, ...)
ov.pl.dynamic_heatmap(adata, *, pseudotime, var_names, ...)
ov.pl.dynamic_trends(res, *, genes, ...)
```

## Source-grounded behavior

**`Monocle` class:**
- 34-line class docstring already documents the workflow + key methods.
- `mono.adata` is a *copy* of the input — modifying upstream `adata` after `Monocle(...)` does not affect the trajectory.
- DDRTree backbone is the canonical Monocle2 algorithm; runs on the high-dimensional ordering-gene subspace.
- `differential_gene_test` runs per-gene VGAM-style tests for non-zero variation along pseudotime; multi-core via `multiprocessing` (or `joblib`); `cores=-1` uses all detected cores.
- BEAM = Branch Expression Analysis Modeling; tests whether a gene's expression as a function of pseudotime differs significantly between branches at the named branch point.

**`dynamic_features`:**
- 80-line Numpy-style docstring already in source — covers the full GAM-fitting API including the multi-AnnData / multi-group dispatch (`data` can be a `dict[name → AnnData]` to fit several datasets at once).
- Default GAM is gaussian-identity; for raw counts use `distribution='gamma', link='log'`.
- `groupby` + `groups` enable branch-aware fitting: each group gets its own fit; `dynamic_trends(compare_groups=True, split_time=..., shared_trunk=True)` then visualises shared trunks.
- `store_raw=True` keeps original cell-level data + the `raw_obs_keys` fields in the result so downstream plotters can colour points.
- Returns a `DynamicFeaturesResult` with `.fitted`, `.confidence`, `.failed` (genes that didn't converge), `.raw` (when `store_raw`).

**`branch_streamplot` / `dynamic_heatmap` / `dynamic_trends`:**
- Plotting layer over the GAM result; signatures verified from notebook usage. The plotting helpers all accept the standard matplotlib `ax`, `figsize`, `show` kwargs.
- `dynamic_heatmap(use_fitted=True)` reads from the GAM-fitted curves; `use_fitted=False` falls back to raw-binned smoothing.
- `dynamic_trends(compare_features=True)` overlays multiple genes on one axis; `compare_groups=True` requires a `groupby` `dynamic_features` result.

## Notebook ↔ skill alignment

| Notebook section | Skill section |
|---|---|
| Olsson WT load | Quick Workflow §1; Input Contract |
| `Monocle().preprocess().select_ordering_genes(max_genes=1000).plot_ordering_genes()` | Quick Workflow §2-3 |
| `mono.plot_trajectory(color_by='State', ...)` | Quick Workflow §5 |
| `ov.pl.branch_streamplot` | Quick Workflow §6 |
| Re-instantiate on ordering subset + `differential_gene_test(cores=-1)` | Quick Workflow §7; reference.md |
| `ov.pl.dynamic_heatmap(top40, use_fitted=True, cell_bins=200, order_by='peak')` | Quick Workflow §8 |
| `ov.single.dynamic_features` + `ov.pl.dynamic_trends` (single-gene + compare_features + compare_groups) | Quick Workflow §9-10 |
| `mono_ord.BEAM(branch_point=1, cores=-1)` + branch-aware visualisation | Quick Workflow §11 |
| `ov.pl.dynamic_heatmap(var_names=panels)` (gene-panel mode) | reference.md (multi-set heatmap block) |

## Docstring supplementation log

`Monocle` (34L), `dynamic_features` (80L) — already well-documented Numpy-style; no supplementation needed.

`differential_gene_test`, `BEAM`, `plot_trajectory`, `select_ordering_genes`, `plot_ordering_genes` — methods on `Monocle`. The class-level docstring documents the workflow; per-method docstrings vary in length but were verified to be present and accurate against the notebook usage.

`ov.pl.branch_streamplot`, `ov.pl.dynamic_heatmap`, `ov.pl.dynamic_trends` — plotting helpers in `ov.pl`. Verified parameter names from notebook calls; full signatures readable in `omicverse/pl/`.

## Reviewer-Run Empirical Checks

- All cited functions importable: `from omicverse.single import Monocle, dynamic_features; from omicverse.pl import branch_streamplot, dynamic_heatmap, dynamic_trends` ✓
- `Monocle.adata` is a copy — verified in source (`__init__` calls `adata.copy()`).
- BEAM hits should be a *subset* of `differential_gene_test` hits with the same q cutoff (constructive: BEAM adds the branch-dependence constraint).
- No live smoke run executed; the Olsson WT dataset (~5k cells) is the canonical smoke target.

## Known caveats

- The `omicverse.single._via.pyVIA` import currently raises `NameError: name 'matplotlib' is not defined` due to a missing local import inside the class body (line ~358 of `_via.py`). Doesn't affect this skill but is flagged in the VIA skill's source-grounding doc.
