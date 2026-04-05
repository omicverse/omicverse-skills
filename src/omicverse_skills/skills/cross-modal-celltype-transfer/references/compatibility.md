# Compatibility Notes

- The notebook uses `X_glue`, but the reusable skill should work with any shared latent space already present on both objects.
- `mde` is optional; it depends on `pymde`, which may be absent in lightweight environments.
- `mode='paper'` is not usable in the current source even though the docstring mentions it.
- `label_keys` is a prefix filter. If you use a broad prefix, the transfer will produce multiple output columns.
- The notebook's `X_umap` plot is not safe to require unless the basis already exists.
- If you want a single output column, store the first prediction column explicitly in `transf_celltype` and its uncertainty in `transf_celltype_unc`.
