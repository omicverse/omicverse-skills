# Compatibility Notes

## Notebook Drift vs Current Wrapper Behavior

- The notebook presents only a subset of the current `cluster(...)` `method` surface.
- Current code writes both `mclust` and `gmm_cluster` for the `GMM` path.
- Current code returns a fitted model only for `scICE`; the other backends mutate `adata.obs` and return `None`.
- The Louvain path still depends on the external Louvain stack and can be more fragile than Leiden in some environments.

## Runtime Considerations

- `scICE` is materially more expensive than Leiden, Louvain, or a small GMM run because it searches over a cluster-count range and uses bootstrap-based consistency scoring.
- Prefer OmicVerse preprocessing wrappers such as `ov.pp.normalize_total`, `ov.pp.log1p`, `ov.pp.pca`, and `ov.pp.neighbors` in reusable examples and smoke utilities.
- Keep acceptance and smoke execution shell-neutral.
