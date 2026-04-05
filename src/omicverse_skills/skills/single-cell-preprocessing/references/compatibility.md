# Compatibility Notes

## Notebook-To-Source Differences

- The notebook uses `ov.settings.cpu_gpu_mixed_init()`, but the reusable skill should treat that as an optional runtime toggle rather than a hard requirement.
- The notebook's `ov.pp.preprocess(..., mode='shiftlog|pearson')` path is preserved, but the live source also exposes `shiftlog|seurat`, `pearson|pearson`, and `pearson|seurat`.
- `ov.pp.pca(adata, layer='scaled', ...)` stores the downstream representation under `scaled|original|X_pca`; use that key in `neighbors` and any plot call that needs the PCA basis.
- `ov.pl.mde(...)` copies `X_mde` to `MDE` by default for plotting. Disable conversion if you need to preserve the original key.
- `ov.single.find_markers(..., method='cosg')` expects raw counts, while statistical methods such as `wilcoxon` expect log-normalized data.

## Smoke Validation Caveat

- Tiny synthetic fixtures can be too small for scrublet after QC filtering. The reviewer smoke path uses a lightweight QC callability check instead of treating scrublet as mandatory.
- The notebook's `recover_counts` cell is diagnostic and is not required for the reusable preprocessing spine.
- The `ConvexHull` and `gen_mpl_labels` cells are visualization polish only; keep them optional unless a caller explicitly wants the same figure styling.

## Current Runtime Notes

- The inspected CPU path is the one validated in the smoke checks.
- The notebook's mixed-mode branch is still valid, but it should be selected only when the runtime really needs it.
