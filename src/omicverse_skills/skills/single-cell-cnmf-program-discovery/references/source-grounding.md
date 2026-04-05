# Source Grounding

## Live Interfaces Checked

- `ov.single.cNMF(...)`
- `cNMF.factorize(self, worker_i=0, total_workers=1)`
- `cNMF.consensus(self, k, density_threshold=0.5, local_neighborhood_size=0.3, show_clustering=True, ...)`
- `cNMF.load_results(self, K, density_threshold, n_top_genes=100, norm_usage=True)`
- `cNMF.get_results(self, adata, result_dict)`
- `cNMF.get_results_rfc(self, adata, result_dict, use_rep='STAGATE', cNMF_threshold=0.5)`

## Source-Verified Branches

- `use_gpu` and `gpu_id` control accelerator selection for factorization.
- `factorize(...)` is worker-aware and requires all workers to run when `total_workers > 1`.
- `load_results(...)` normalizes usages when `norm_usage=True` and renames columns to `cNMF_<k>`.
- `get_results(...)` writes `cNMF_cluster` from max usage.
- `get_results_rfc(...)` writes `cNMF_cluster_rfc` and `cNMF_cluster_clf` from tree-based classifiers trained on thresholded usage subsets.

## Practical Branch Guidance

- Use single-worker CPU mode for bounded smoke and local review.
- Use multi-worker mode only when the caller can orchestrate every worker index.
- Use the RFC path only after normalized usages have already been loaded.

## Reviewer-Side Evidence

- The notebook was partitioned into graph clustering, LDA, cNMF, and evaluation before writing the skill.
- Live signatures and source behavior were checked for constructor, factorization, consensus, load, and labeling methods.
