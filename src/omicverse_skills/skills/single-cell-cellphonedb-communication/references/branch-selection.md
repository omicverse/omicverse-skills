# Branch Selection

## Notebook Capability Partition

- Core executable job: run CellPhoneDB statistical analysis and produce a visualization-ready interaction `AnnData`.
- Optional downstream jobs:
  - aggregated sender-receiver communication network
  - pathway-level communication scoring
  - ligand-receptor contribution analysis
  - bubble and chord visualizations
  - signaling-role centrality and heatmaps
- Notebook-only pedagogy: CellPhoneDB background and interpretation text.

This stays as one skill because all downstream jobs consume the same `adata_cpdb` object and `CellChatViz` state. Splitting by plot type would create thin wrappers with no stable independent input contract beyond that object.

## `run_cellphonedb_v5(...)` Branches

### Cell filtering

- `min_cell_fraction` filters cell types before analysis.
- `min_genes` and `min_cells` are forwarded into the wrapper’s `scanpy` filtering stage.

### Statistical depth

- `iterations` controls permutation count.
- `pvalue` controls the significance threshold used in CellPhoneDB result generation.
- `threshold` controls the minimum fraction of cells expressing a gene for analysis.

### Runtime and file-management behavior

- `threads` sets CellPhoneDB thread count.
- `output_dir=None` creates an output directory automatically; supplying a directory keeps results in a known place.
- `temp_dir=None` creates a temporary working directory automatically.
- `cleanup_temp=True` removes wrapper-created temporary input files after the run.
- `debug=True` requests CellPhoneDB intermediate tables.
- `separator` controls how sender-receiver pair columns are parsed into the visualization `AnnData`.

### Additional CellPhoneDB branch surface

- `**kwargs` are forwarded into the underlying CellPhoneDB statistical call, so this wrapper can expose more CellPhoneDB options than the notebook uses.

## `CellChatViz` and Derived Branches

### Pathway aggregation

`compute_pathway_communication(method=...)` supports:

- `mean`
- `sum`
- `max`
- `median`

### Aggregate layout

`netVisual_aggregate(layout=...)` supports:

- `circle`
- `hierarchy`

### Chord view focus

`netVisual_chord_cell(...)` branches on:

- `signaling=None` versus one or more specific pathways
- `group_celltype=None` versus grouped cell-type mapping
- `sources` and `targets` filters
- `normalize_to_sender=True` versus `False`

### Bubble view encoding

`netVisual_bubble_marsilea(...)` branches on:

- `signaling=None` versus pathway-specific views
- `sources_use` and `targets_use`
- `group_pathways=True` versus `False`
- `transpose=True` versus `False`
- `scale=None`, `row`, `column`, `row_minmax`, or `column_minmax`
- sender and receiver color-bar toggles

### Signaling-role views

`netAnalysis_signalingRole_heatmap(pattern=...)` supports:

- `outgoing`
- `incoming`

It also branches on `row_scale` and optional pathway selection.

## Practical Default Path

- Use CellPhoneDB analysis first.
- Move next to `compute_aggregated_network(...)` for a stable high-level summary.
- Use pathway aggregation only when the user asks for pathway-level prioritization.
- Use bubble or chord plots only when the user wants a focused subset of pathways or ligand-receptor pairs.
- Use signaling-role centrality and heatmaps when the user asks who acts as sender, receiver, mediator, or dominant signaling hub.
