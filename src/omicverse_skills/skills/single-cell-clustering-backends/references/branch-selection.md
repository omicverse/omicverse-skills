# Branch Selection

## Backend Choice

- Choose `leiden` as the default graph clustering path.
- Choose `louvain` when Louvain compatibility matters more than using the newer default.
- Choose `GMM` when `n_components` is known and the user wants embedding-based Gaussian mixture clustering.
- Choose `scICE` when the user wants candidate stable solutions across a cluster-count range.

## Backend-Specific Inputs

- `leiden` and `louvain`: require an existing graph.
- `GMM`: requires `use_rep` and `n_components`; `covariance_type` materially changes the fitted model.
- `scICE`: requires `use_rep` and tuning of `resolution_range`, `n_boot`, and `n_steps`.

## Output Interpretation

- `GMM` produces both `mclust` and `gmm_cluster`.
- `scICE` may produce multiple `scICE_k*` outputs rather than one canonical cluster key.
