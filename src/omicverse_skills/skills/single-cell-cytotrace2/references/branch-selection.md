# Branch Selection

## Notebook Capability Partition

- Core executable job: preprocess expression values if needed, stage pretrained weights, run `ov.single.cytotrace2(...)`, and capture potency outputs.
- Optional downstream job: overlay the potency outputs on an existing embedding.
- Notebook-only pedagogy: biological explanation of potency categories and visual interpretation.

One skill is the better boundary here because the plotting stage is only a thin consumer of the CytoTRACE2 output columns and does not have a distinct input contract.

## `ov.pp.preprocess(...)` Branches

The live implementation splits `mode` on `|` and evaluates each token separately.

### Normalization token

- `shiftlog`: `normalize_total(...)` then `log1p(...)`.
- `pearson`: `normalize_pearson_residuals(...)`.

### HVG token

- `pearson`: `highly_variable_genes(..., flavor='pearson_residuals', layer='counts')`.
- `seurat`: `highly_variable_genes(..., flavor='seurat_v3', layer='counts')`.

Additional preprocessing knobs that materially affect behavior:

- `organism`: matters only when `no_cc=True` triggers cell-cycle-gene removal.
- `batch_key`: passed into HVG selection.
- `identify_robust`: controls the up-front robust-gene filter and can change the gene universe before HVG selection.

## `ov.single.cytotrace2(...)` Branches

### Species branch

- `mouse`
- `human`

The wrapper does not hard-enforce the species choice, but it does warn when gene-symbol casing suggests the other species.

### Parallelization branch

- `disable_parallelization=False`: worker counts are derived from chunk counts, smoothing chunk counts, detected CPU cores, and optional `max_cores`.
- `disable_parallelization=True`: forces one prediction worker and one smoothing worker.

### Batching branch

- `batch_size`: controls prediction chunking.
- `smooth_batch_size`: controls smoothing chunking inside each prediction chunk.
- `max_pcs`: passed to KNN smoothing.

### Small-dataset edge branch

- When a chunk has fewer than 100 cells, `process_subset(...)` skips KNN smoothing and returns the pre-smoothing table.
- The outer wrapper later expects `CytoTRACE2_Score`, which is missing on that branch, so current OmicVerse raises `KeyError: 'CytoTRACE2_Score'`.
- For a working smoke path, keep the processed chunk size at 100 cells or more.

## Visualization Branch

- The notebook uses `ov.utils.embedding(...)`.
- Current OmicVerse exposes that as a compatibility wrapper, so treat it as a presentational branch and validate the potency columns first.
