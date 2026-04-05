# Source Grounding

## Interfaces Checked

The skill content was grounded against live OmicVerse code in the `omicverse` conda environment using:

- `inspect.signature(...)`
- `inspect.getdoc(...)`
- `inspect.getsource(...)`
- a reviewer-run bounded CytoTRACE2 smoke execution

## `ov.pp.preprocess`

Observed signature:

```python
(adata, mode='shiftlog|pearson', target_sum=500000.0, n_HVGs=2000, organism='human', no_cc=False, batch_key=None, identify_robust=True)
```

Important source-grounded behavior:

- Saves the original matrix to `adata.layers['counts']`.
- Applies a robust-gene filter before normalization and HVG selection when `identify_robust=True`.
- Splits `mode` on `|` and branches separately for normalization and HVG selection.
- Supports `shiftlog` and `pearson` in the normalization slot.
- Supports `pearson` and `seurat` in the HVG slot.
- Mirrors the processed state back to the original object if an internal slice was created.

## `ov.single.cytotrace2`

Observed signature:

```python
(adata, use_model_dir, species='mouse', batch_size=10000, smooth_batch_size=1000, disable_parallelization=False, max_cores=None, max_pcs=200, seed=14, output_dir='cytotrace2_results')
```

Important source-grounded behavior:

- Reads expression values from `adata.X`.
- Converts sparse matrices with `.toarray()` before building the prediction frame.
- Warns when gene-symbol casing looks inconsistent with the declared species.
- Resets `batch_size` to the dataset size when the requested batch is larger than the cell count.
- Writes the final five CytoTRACE2 output columns into `adata.obs`.
- Writes a tab-delimited result table to the output directory.

## Parallelization Helper

`calculate_cores_to_use(...)` was inspected directly.

Observed behavior:

- If both chunk counts are 1, it prints that prediction will not be parallelized.
- If parallelization is enabled and `max_cores` is unset, worker counts are derived from detected CPUs.
- If `max_cores` is set, it caps prediction and smoothing workers from that limit.

## `process_subset(...)`

Important source-grounded behavior:

- Calls CytoTRACE2 preprocessing, top-variable-gene selection, prediction, and smoothing helpers from the bundled external module.
- Writes intermediate ranked, binned, and top-variable-gene files per chunk.
- Skips KNN smoothing entirely when the chunk has fewer than 100 cells.

## Plotting Surface

`ov.utils.embedding` was inspected and currently resolves through a wrapper signature of `(*args, **kwargs)`. The notebook plotting call is therefore documented as an optional overlay stage, not as a strongly typed core interface.

## Reviewer-Run Empirical Checks

- A 60-cell local smoke run reproduced the current small-dataset failure path: `KeyError: 'CytoTRACE2_Score'`.
- A 120-cell local smoke run completed successfully on local dentate gyrus data after `ov.pp.preprocess(...)`.
- The successful run wrote the expected five CytoTRACE2 columns into `adata.obs` and produced the final result table.
