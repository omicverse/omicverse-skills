# Source Grounding

This skill is grounded in the local OmicVerse source tree and in live interface inspection under the `omictest` conda environment.

## Primary Source Files

- `omicverse/pp/_qc.py`
- `omicverse/pp/_preprocess.py`
- `omicverse/pp/_pca.py`
- `omicverse/pp/_neighbors.py`
- `omicverse/pp/_umap.py`
- `omicverse/pp/_leiden.py`
- `omicverse/single/_markers.py`
- `omicverse/pl/_dotplot.py`

## Inspected Interfaces

### `ov.pp.preprocess`

Signature:

```text
(adata, mode='shiftlog|pearson', target_sum=500000.0, n_HVGs=2000, organism='human', no_cc=False, batch_key=None, identify_robust=True)
```

Observed source branches:

- `mode='shiftlog|pearson'`
- `mode='shiftlog|seurat'`
- `mode='pearson|pearson'`
- `mode='pearson|seurat'`

### `ov.pp.qc`

The public wrapper is `qc(adata, **kwargs)`, and the live source dispatches through CPU, mixed, or GPU paths.

Observed QC options in source:

- `mode='seurat'`
- `mode='mads'`
- `doublets_method='scrublet'`
- `doublets_method='sccomposite'`
- `filter_doublets=True|False`

### `ov.pp.pca`

Signature:

```text
(adata, n_pcs=50, layer='scaled', inplace=True, **kwargs)
```

The notebook's `layer='scaled'` path stores PCA under:

- `adata.obsm['scaled|original|X_pca']`
- `adata.varm['scaled|original|pca_loadings']`
- `adata.uns['scaled|original|pca_var_ratios']`

### `ov.pp.neighbors`

Signature:

```text
(adata: AnnData, n_neighbors: int = 15, n_pcs: Optional[int] = None, use_rep: Optional[str] = None, knn: bool = True, random_state: int = 0, n_jobs: Optional[int] = None, method: Optional[Literal['umap', 'gauss', 'rapids']] = 'umap', transformer: Optional[str] = None, metric: Union[... ] = 'euclidean', metric_kwds: Mapping[str, Any] = mappingproxy({}), key_added: Optional[str] = None, copy: bool = False, **kwargs)
```

Observed branch values:

- `method='umap'`
- `method='gauss'`
- `method='rapids'`
- `transformer='pyg'`
- `transformer='pynndescent'`
- `transformer='rapids'`

### `ov.pp.leiden`

Signature:

```text
(adata, resolution=1.0, random_state=0, key_added='leiden', local_iterations=100, max_levels=10, device='cpu', symmetrize=None, **kwargs)
```

Source behavior:

- CPU mode routes to Scanpy-backed Leiden.
- `cpu-gpu-mixed` routes to the omicverse mixed implementation.
- GPU mode routes to RAPIDS where available.

### `ov.single.find_markers`

Signature:

```text
(adata: 'AnnData', groupby: 'str', method: 'str' = 'cosg', n_genes: 'int' = 50, key_added: 'Optional[str]' = None, use_raw: 'Optional[bool]' = None, layer: 'Optional[str]' = None, groups: 'Union[str, Sequence[str]]' = 'all', reference: 'str' = 'rest', corr_method: 'str' = 'benjamini-hochberg', rankby_abs: 'bool' = False, tie_correct: 'bool' = False, pts: 'bool' = True, **kwargs) -> 'None'
```

Observed branch values:

- `method='cosg'`
- `method='wilcoxon'`
- `method='t-test'`
- `method='t-test_overestim_var'`
- `method='logreg'`

Source behavior:

- `cosg` is the fast raw-count branch.
- statistical methods use log-normalized data.

### `ov.single.get_markers`

Signature:

```text
(adata: 'AnnData', n_genes: 'int' = 10, key: 'str' = 'rank_genes_groups', groups: 'Optional[Union[str, Sequence[str]]]' = None, return_type: 'str' = 'dataframe', min_logfoldchange: 'Optional[float]' = None, min_score: 'Optional[float]' = None, min_pval_adj: 'Optional[float]' = None) -> 'Union[pd.DataFrame, dict]'
```

Observed branch values:

- `return_type='dataframe'`
- `return_type='dict'`

### `ov.pl.markers_dotplot`

Signature:

```text
(adata: anndata.AnnData, groupby: Optional[str] = None, key: Optional[str] = None, n_genes: int = 5, groups: Union[str, Sequence[str], NoneType] = None, standard_scale: Optional[Literal['var', 'group']] = 'var', cmap: Union[matplotlib.colors.Colormap, str, NoneType] = 'Spectral_r', dendrogram: bool = False, min_logfoldchange: Optional[float] = None, use_raw: Optional[bool] = None, layer: Optional[str] = None, figsize: Optional[Tuple[float, float]] = None, show: Optional[bool] = None, save: Union[bool, str, NoneType] = None, return_fig: bool = False, **kwds) -> Optional[Any]
```

Source behavior:

- `show=False` is the useful smoke path when you only want to verify the plotting branch without opening a window.

## Interface Inspection Evidence

- `inspect.signature(...)` was run against the live OmicVerse checkout through the `omictest` Python runtime.
- The notebook's branch-heavy parameters were checked against source, not only the rendered tutorial cells.
