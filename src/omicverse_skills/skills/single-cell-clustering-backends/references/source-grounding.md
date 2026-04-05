# Source Grounding

## Live Interfaces Checked

- `ov.utils.cluster(adata, method='leiden', use_rep='X_pca', random_state=1024, n_components=None, **kwargs)`

## Source-Verified Branches

- Current code supports `leiden`, `louvain`, `kmeans`, `GMM`, `mclust`, `mclust_R`, `schist`, and `scICE`.
- The notebook only demonstrates `leiden`, `louvain`, `scICE`, and `GMM`; it does not exhaust the current API.
- `leiden` and `louvain` forward to scanpy graph clustering.
- OmicVerse exposes `ov.pp.normalize_total(...)`, `ov.pp.log1p(...)`, `ov.pp.pca(...)`, and `ov.pp.neighbors(...)`; reusable examples in this skill should prefer those wrappers instead of dropping to raw Scanpy preprocessing calls.
- `GMM` and `mclust` call the same Gaussian-mixture helper and write both `mclust` and `gmm_cluster`.
- `scICE` instantiates a `scICE` model with `use_gpu=False` inside this wrapper, returns that model, and writes `scICE_k*` columns through `add_to_adata(...)`.
- `scICE` stable-label columns depend on the discovered cluster numbers and are not fixed to one suffix such as `scICE_k13`.

## Practical Branch Guidance

- Use graph-based clustering when the graph is already part of the analysis contract.
- Use `GMM` when the user wants mixture-model clustering on a chosen embedding.
- Use `scICE` when the user wants consistency-based cluster-number exploration rather than one single graph partition.

## Reviewer-Side Evidence

- The notebook was partitioned into graph clustering, topic modeling, cNMF, and evaluation before writing the skill.
- The installed OmicVerse source was inspected for current `method` branches and backend-specific output columns.
