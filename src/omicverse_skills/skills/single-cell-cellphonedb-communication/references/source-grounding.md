# Source Grounding

## Interfaces Checked

The skill content was grounded against live OmicVerse and CellPhoneDB code in the `omicverse` conda environment using:

- `inspect.signature(...)`
- `inspect.getdoc(...)`
- `inspect.getsource(...)`
- the repository interface-inspection helper
- a reviewer-run bounded CellPhoneDB smoke path on local hematopoiesis data

## `ov.single.run_cellphonedb_v5`

Observed signature:

```python
(adata, cpdb_file_path, celltype_key='celltype', min_cell_fraction=0.005, min_genes=200, min_cells=3, iterations=1000, threshold=0.1, pvalue=0.05, threads=10, output_dir=None, temp_dir=None, cleanup_temp=True, debug=False, separator='|', **kwargs)
```

Important source-grounded behavior:

- validates or auto-downloads the CellPhoneDB database archive
- checks that `celltype_key` exists in `adata.obs`
- filters rare cell types by `min_cell_fraction`
- uses `adata.raw` if present, otherwise `adata.X`
- runs `scanpy` cell/gene filtering before CellPhoneDB
- writes temporary counts and metadata files for CellPhoneDB
- forwards extra keyword arguments into `cpdb_statistical_analysis_method.call(...)`
- converts CellPhoneDB results into a visualization-ready `AnnData` through `format_cpdb_results_for_viz(...)`

## `format_cpdb_results_for_viz(...)`

Observed signature:

```python
(cpdb_results, separator='|')
```

Important source-grounded behavior:

- identifies sender-receiver pair columns by the configured separator
- transposes pair columns so cell-type pairs become observations
- stores `means` and `pvalues` in layers
- writes `sender` and `receiver` into `obs`
- preserves interaction metadata in `var`

## `CellChatViz` Core

Observed constructor signature:

```python
(adata, palette=None)
```

Important source-grounded behavior:

- accepts either explicit palette dicts, palette sequences, or no palette
- infers unique cell types from `sender` and `receiver`
- fills missing colors from existing metadata or fallback palettes

## Branch-Heavy Downstream Methods

### `compute_pathway_communication`

Observed signature:

```python
(self, method='mean', min_lr_pairs=1, min_expression=0.1)
```

Observed branch values:

- `mean`
- `sum`
- `max`
- `median`

### `netVisual_aggregate`

Observed signature:

```python
(self, signaling, layout='circle', vertex_receiver=None, vertex_sender=None, pvalue_threshold=0.05, ...)
```

Observed branch values:

- `circle`
- `hierarchy`

### `netVisual_bubble_marsilea`

Observed signature includes pathway selection, sender/receiver filtering, grouping, transposition, and scaling.

Observed scaling branch values:

- `None`
- `row`
- `column`
- `row_minmax`
- `column_minmax`

### `netAnalysis_signalingRole_heatmap`

Observed signature:

```python
(self, pattern='outgoing', signaling=None, row_scale=True, figsize=(12, 8), cmap='RdYlBu_r', show_totals=True, title=None, save=None, min_threshold=0.1)
```

Observed pattern branch values:

- `outgoing`
- `incoming`

## Reviewer-Run Empirical Checks

- A reduced local run with `iterations=10` failed inside CellPhoneDB with `ZeroDivisionError`, so very small permutation counts are not safe as a smoke default.
- A bounded local run with `iterations=100`, `threads=1`, and a three-cell-type subset of local hematopoiesis data completed successfully.
- The successful run produced:
  - a non-empty CellPhoneDB result dict
  - a visualization-ready `adata_cpdb`
  - `means` and `pvalues` layers
  - `sender` and `receiver` observation metadata
  - working aggregated-network, pathway-communication, and centrality calculations
