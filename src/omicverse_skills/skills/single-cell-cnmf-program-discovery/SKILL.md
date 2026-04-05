---
name: omicverse-single-cell-cnmf-program-discovery
description: Run OmicVerse single-cell cNMF program discovery as a reusable, triggerable skill. Use when fitting consensus NMF gene programs on single-cell AnnData, choosing between CPU and GPU factorization, or converting normalized usage programs into hard cluster labels and RFC-based labels.
---

# OmicVerse Single-Cell cNMF Program Discovery

## Goal

Turn the notebook's cNMF section into a reusable job: run consensus NMF over candidate ranks, choose one rank, load normalized program usages, and optionally derive hard labels or RFC-based labels.

## Quick Workflow

1. Inspect whether the data matrix and PCA embedding are ready for downstream labeling.
2. Choose candidate `components`, `n_iter`, whether to use GPU, and whether results should be persisted to an output directory.
3. Run `factorize(...)` and `combine(...)`.
4. Run `consensus(...)`, then `load_results(...)`.
5. Use `get_results(...)` for direct max-usage labels or `get_results_rfc(...)` for classifier-derived labels.

## Interface Summary

- `ov.single.cNMF(adata, components, n_iter=100, densify=False, tpm_fn=None, seed=None, beta_loss='frobenius', num_highvar_genes=2000, genes_file=None, alpha_usage=0.0, alpha_spectra=0.0, init='random', output_dir=None, name=None, use_gpu=True, gpu_id=0)` constructs the workflow wrapper.
- `factorize(worker_i=0, total_workers=1)` runs NMF iterations for the worker's assigned jobs.
- `combine(skip_missing_files=False)` merges replicate factorizations.
- `consensus(k, density_threshold=0.5, local_neighborhood_size=0.3, show_clustering=True, ...)` selects a rank and produces consensus program usage.
- `load_results(K, density_threshold, n_top_genes=100, norm_usage=True)` returns normalized usage and top-gene summaries.
- `get_results(adata, result_dict)` writes `cNMF_cluster`.
- `get_results_rfc(adata, result_dict, use_rep='STAGATE', cNMF_threshold=0.5)` writes `cNMF_cluster_rfc` and `cNMF_cluster_clf`.

## Stage Selection

- Use CPU factorization for bounded smoke runs or when GPU is unavailable.
- Use GPU only when the environment really supports it and the workload justifies it.
- Use `get_results(...)` for direct max-usage labeling.
- Use the RFC path only when the user explicitly wants classifier-derived hard labels on top of the usage matrix.
- Treat K-selection figures and embedding panels as optional reporting, not the core skill contract.

## Input Contract

- Start from `AnnData`.
- Choose candidate ranks in `components`.
- Provide a writable output location if you need persisted intermediate files.
- Ensure the embedding named by `use_rep` exists before the RFC path.

## Minimal Execution Patterns

```python
import numpy as np
import omicverse as ov

cnmf = ov.single.cNMF(
    adata,
    components=np.arange(5, 11),
    n_iter=20,
    seed=14,
    num_highvar_genes=2000,
    output_dir="...",
    name="dg_cNMF",
    use_gpu=False,
)
cnmf.factorize(worker_i=0, total_workers=1)
cnmf.combine(skip_missing_files=True)
```

```python
cnmf.consensus(k=7, density_threshold=2.0, show_clustering=True)
result_dict = cnmf.load_results(K=7, density_threshold=2.0)
cnmf.get_results(adata, result_dict)
```

```python
cnmf.get_results_rfc(
    adata,
    result_dict,
    use_rep="scaled|original|X_pca",
    cNMF_threshold=0.5,
)
```

## Constraints

- Do not pretend multi-worker factorization is complete unless every worker ran.
- Use `use_gpu=True` only when accelerator support is actually available and intended.
- Do not assume the notebook's chosen `k` or `density_threshold` transfers to other datasets.
- Do not run the RFC path without a real embedding in `obsm`.
- Keep smoke and acceptance commands shell-agnostic.

## Validation

- Check that `factorize(...)` actually ran the assigned iterations.
- Check that `consensus(...)` completed before `load_results(...)`.
- After `load_results(...)`, check that normalized usage columns exist.
- After `get_results(...)`, check that `cNMF_cluster` exists.
- After the RFC path, check that `cNMF_cluster_rfc` and `cNMF_cluster_clf` exist.
- If only a bounded smoke path was run, say which expensive stages were reduced.

## Resource Map

- Use the branch selection notes when choosing CPU vs GPU, single-worker vs multi-worker, or direct labels vs RFC labels.
- Use the source grounding notes for current signatures, worker semantics, and output labels.
- Use the notebook mapping notes to trace the notebook's cNMF section into this reusable skill.
- Use the compatibility notes for compute and filesystem-sensitive behavior.
