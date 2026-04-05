# Source Notebook Map

## Notebook sections to reusable steps

- Load RNA reference and ATAC query: create the paired input contract.
- Optional concatenation: only used for joint inspection, not for transfer itself.
- Optional `mde` embedding: display branch only.
- Weighted KNN trainer: core model fit on the reference embedding.
- Weighted KNN transfer: core label transfer onto the query.
- Writeback to `obs`: store transfer labels and uncertainty.
- Final plots: confirm the transferred labels on any existing basis.

## Notebook-only details

- The tutorial's file paths are sample data, not part of the reusable skill.
- The exact `X_mde` plots are presentation choices, not required outputs.
- The final `X_umap` plot is only valid if that basis already exists.
