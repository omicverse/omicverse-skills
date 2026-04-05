# Notebook Mapping Notes

## Reused upstream prerequisite

The notebook's QC and preprocessing cells are treated as upstream prerequisites, not as the core identity of this skill:

- QC with `ov.pp.qc(...)`
- preprocessing with `ov.pp.preprocess(...)`
- scaling and PCA with `ov.pp.scale(...)` and `ov.pp.pca(...)`

Those belong to the existing preprocessing skill and are intentionally not expanded into a second preprocessing skill here.

## Stage 1: batch-labeled merged input

Notebook cells:

- concatenate multiple samples
- create `adata.obs['batch']`
- preserve integer-like counts before downstream integration

Reusable contract extracted from that stage:

- an `AnnData` object with batch labels already present

## Stage 2: integration backends

Notebook cells map directly to `methods=` branches:

- `methods='harmony'`
- `methods='combat'`
- `methods='scanorama'`
- `methods='scVI'`
- `methods='CellANOVA'`
- `methods='Concord'`

These remain inside one skill because they are alternative branches of one callable with the same high-level job: batch integration.

## Stage 3: optional embedding visualization

Notebook pattern:

- compute an MDE projection from an integrated embedding
- color by `batch` and `cell_type`

This is kept as an optional stage because it only becomes meaningful after an integration embedding exists.

## Stage 4: optional benchmarking

Notebook pattern:

- instantiate `Benchmarker(...)`
- supply a non-integrated reference embedding plus multiple integrated embeddings
- run `.benchmark()`
- optionally call `.plot_results_table(...)`

This stage stays inside the same skill because it is a direct downstream evaluation of the integration outputs and is not usually useful without them.
