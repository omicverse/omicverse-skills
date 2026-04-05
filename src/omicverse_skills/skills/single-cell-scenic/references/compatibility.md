# Compatibility

## Current Runtime Notes

- The live `SCENIC.cal_grn(...)` default layer is `counts`, while the notebook example calls `layer='raw_count'`.
- The GRN step expects a raw-count sparse layer and converts it internally with `.toarray()`.
- Full regulon pruning depends on external cisTarget ranking databases and a motif annotation table; the skill does not assume those files can be downloaded on demand.
- `rho_mask_dropouts=True` recreates the notebook-era dropout masking behavior, while the installed pySCENIC default now uses all cells.
- The current wrapper's arboreto-backed `grnboost2` path failed on a synthetic reviewer smoke because the backend received expression values without gene names. Validate `grnboost2` and `genie3` explicitly before relying on them.
- The optional TF-network plotting stage also depends on `adjustText`, `networkx`, and `umap`.

## Validation Scope

- A lightweight reviewer smoke can validate constructor behavior plus a bounded `regdiffusion` branch without claiming a full cisTarget reproduction.
- A full end-to-end regulon run remains more expensive because it combines module generation, motif pruning, and AUCell scoring against external resources.
- If the user wants exact tutorial figures, keep the worked-example grouping column and TF list explicit instead of assuming they generalize.
