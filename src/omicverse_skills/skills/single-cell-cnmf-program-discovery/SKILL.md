---
name: omicverse-single-cell-cnmf-program-discovery
description: Run OmicVerse single-cell NMF program discovery as a reusable, triggerable skill — both the classical Python `ov.single.cNMF` (consensus NMF with CPU/GPU factorization, K-selection, RFC labelling) and the Rust-backed `ov.single.NMF` (fast `nmf-rs` backend: dnmf default, Brunet-style K-selection with stability-drop auto-K, cNMF-style consensus heatmap, RFC labels). Use when fitting consensus NMF gene programs on single-cell AnnData, choosing K, building consensus, or converting normalized usage programs into hard cluster labels.
---

# OmicVerse Single-Cell NMF / cNMF Program Discovery

## Goal

Turn the notebook's NMF section into a reusable job. Two backends are supported and share the same downstream contract (normalized usage matrix + per-program top genes + optional RFC labels):

- **`ov.single.cNMF`** — classical Kotliar-style consensus NMF; CPU or `torchnmf` GPU; multi-worker `factorize → combine → consensus → load_results`. Use when reproducing the cNMF paper pipeline exactly or when you need its K-selection / silhouette diagnostics.
- **`ov.single.NMF`** — Rust port via the optional `nmf-rs` package (`pip install nmf-rs`); ~75–280× faster than R-`NMF`. Default recipe is `method='dnmf'` + NNDSVD init + 25 iterations. Provides Brunet-style stability + reconstruction K-selection with an auto-detected K via a stability-drop heuristic (Brunet 2004 + Kim-Park 2007), cNMF-style consensus heatmap (`n_runs=50`), and the same RFC label path. Use for fast exploratory NMF, large atlases (>500k cells), or when the cNMF compute envelope is overkill.

## Quick Workflow

**Backend A — `ov.single.cNMF` (classical, multi-worker)**

1. Inspect whether the data matrix and PCA embedding are ready for downstream labeling.
2. Choose candidate `components`, `n_iter`, whether to use GPU, and whether results should be persisted to an output directory.
3. Run `factorize(...)` and `combine(...)`.
4. Run `consensus(...)`, then `load_results(...)`.
5. Use `get_results(...)` for direct max-usage labels or `get_results_rfc(...)` for classifier-derived labels.

**Backend B — `ov.single.NMF` (Rust `nmf-rs`, single-process, fast)**

1. Confirm `nmf-rs` is installed (`pip install nmf-rs`); the wrapper imports it lazily.
2. Instantiate with a placeholder `rank` and run `select_k_brunet(np.arange(5, 16), method='dnmf', n_runs=30, max_iter=50)`; read `auto_k` for the stability-drop choice and inspect with `k_selection_plot(ax=...)`.
3. Re-instantiate `ov.single.NMF(adata, rank=auto_k, ...)` and call `fit(method='dnmf', init='nndsvd', max_iter=25)`.
4. Run `consensus(n_runs=50, method='dnmf', max_iter=25)` and visualise with `plot_consensus_heatmap(...)`.
5. Call `get_results(adata, key_added='NMF', n_top_genes=30)` to write `obsm['NMF_usage']`, `varm['NMF_genes']`, `obs['NMF_module']` and return the usual `{usage_norm, gep_scores, top_genes}` dict.
6. Optional: `get_results_rfc(adata, result_dict, use_rep='scaled|original|X_pca', threshold=0.3, key_added='NMF_module_rfc')` and `plot_top_genes(n_top=10, ...)`.

## Interface Summary

**`ov.single.cNMF` (Backend A)**

- `ov.single.cNMF(adata, components, n_iter=100, densify=False, tpm_fn=None, seed=None, beta_loss='frobenius', num_highvar_genes=2000, genes_file=None, alpha_usage=0.0, alpha_spectra=0.0, init='random', output_dir=None, name=None, use_gpu=True, gpu_id=0)` constructs the workflow wrapper.
- `factorize(worker_i=0, total_workers=1)` runs NMF iterations for the worker's assigned jobs.
- `combine(skip_missing_files=False)` merges replicate factorizations.
- `consensus(k, density_threshold=0.5, local_neighborhood_size=0.3, show_clustering=True, ...)` selects a rank and produces consensus program usage.
- `k_selection_plot(close_fig=False)`, `calculate_silhouette_k(k, density_threshold)`, `plot_silhouette_for_k(...)` — K diagnostics.
- `load_results(K, density_threshold, n_top_genes=100, norm_usage=True)` returns normalized usage and top-gene summaries.
- `get_results(adata, result_dict)` writes `cNMF_cluster`.
- `get_results_rfc(adata, result_dict, use_rep='STAGATE', cNMF_threshold=0.5)` writes `cNMF_cluster_rfc` and `cNMF_cluster_clf`.

**`ov.single.NMF` (Backend B — Rust)**

- `ov.single.NMF(adata, rank, use_hvg=True, num_threads=None, ...)` constructs the fast NMF wrapper. `rank` is a placeholder when you're about to call `select_k_brunet`.
- `.select_k_brunet(k_range, method='dnmf', n_runs=30, max_iter=50, ...) → pd.DataFrame` runs each candidate K many times, returns per-K silhouette + reconstruction; populates `.auto_k` from the stability-drop heuristic (Brunet 2004 + Kim-Park 2007 local-peak rule).
- `.auto_k` — chosen K after `select_k_brunet`.
- `.k_selection_plot(ax=None)` — silhouette (left axis) + reconstruction loss (right axis), with a dashed line at `auto_k`.
- `.fit(method='dnmf', init='nndsvd', max_iter=25)` — main factorisation. `method` ∈ {`'dnmf'` (RcppML-style 2024, default), `'lee'`, `'brunet'`, `'snmf/r'`, `'snmf/l'`, `'ls-nmf'`}; `init` ∈ {`'nndsvd'` (default), `'random'`, ...}.
- `.consensus(n_runs=50, method='dnmf', max_iter=25, ...)` — cNMF-style Brunet consensus over multiple random inits.
- `.plot_consensus_heatmap(figsize=(6, 5))` — averaged binary co-cluster matrix re-ordered by hierarchical clustering of `1 − C̄`.
- `.get_results(adata, key_added='NMF', n_top_genes=30) → dict` writes `obsm[f'{key_added}_usage']`, `varm[f'{key_added}_genes']`, `obs[f'{key_added}_module']`; returns `{usage_norm, gep_scores, top_genes}` (same keys as cNMF `result_dict`).
- `.plot_top_genes(n_top=10, figsize=(8, 7))` — heatmap of top genes per factor.
- `.get_results_rfc(adata, result_dict, use_rep='scaled|original|X_pca', threshold=0.3, key_added='NMF_module_rfc')` — same RFC contract as `cNMF.get_results_rfc` but with a unified `key_added`.

## Stage Selection

**Choose the backend first.**
- Use **`cNMF`** when (a) you need bit-equivalence with the Kotliar cNMF paper, (b) you want multi-worker factorisation across nodes, (c) you specifically need its silhouette K-selection plots.
- Use **`NMF`** (Rust) when (a) you want fast exploratory NMF, (b) you're on a large atlas (>500k cells), (c) you want the `auto_k` stability-drop heuristic, (d) you want a drop-in replacement that yields the same `usage_norm` / `top_genes` / RFC labels.
- Both backends produce comparable usage/top-genes outputs; downstream embedding plots and dot plots use the same `ov.pl.embedding` / `ov.pl.dotplot` calls regardless.

**cNMF-specific:**
- Use CPU factorization for bounded smoke runs or when GPU is unavailable.
- Use GPU only when the environment really supports it (`torchnmf` installed, accelerator available).
- Use `get_results(...)` for direct max-usage labeling.
- Use the RFC path only when the user explicitly wants classifier-derived hard labels on top of the usage matrix.

**NMF (Rust)-specific:**
- Default recipe: `method='dnmf'` + `init='nndsvd'` + `max_iter=25`. Strong biology + fastest.
- Use `method='lee'` / `'brunet'` / `'snmf/r'` / `'snmf/l'` when bit-equivalence with R's `NMF::nmf(...)` is required.
- Use `method='ls-nmf'` with `weight=mask` when missing-value imputation is the actual goal.
- Trust `auto_k` only when the silhouette curve shows a clear plateau-then-drop; on monotonic curves, bump `n_runs` and re-run, or fall back to manual inspection of `k_selection_plot`.
- `n_runs=50` is the canonical setting for `consensus(...)` — enough to see K-block structure on most cohorts.

## Input Contract

- Start from `AnnData` with log-normalised counts in `.X` (NMF requires non-negative input — never pass scaled or PCA data).
- For `cNMF`: choose candidate ranks in `components`; provide a writable `output_dir` if you need persisted intermediate files; ensure the embedding named by `use_rep` exists before the RFC path.
- For `NMF` (Rust): `nmf-rs` installed; ensure the embedding named by `use_rep` exists before `get_results_rfc`.

## Minimal Execution Patterns

**cNMF (classical):**

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

**NMF (Rust `nmf-rs`):**

```python
import numpy as np
import matplotlib.pyplot as plt
import omicverse as ov

# 1) K-selection with Brunet stability + reconstruction
rust_nmf = ov.single.NMF(adata, rank=12, use_hvg=True, num_threads=8)
k_df = rust_nmf.select_k_brunet(
    np.arange(5, 16),
    method='dnmf', n_runs=30, max_iter=50,
)
print('auto_k =', rust_nmf.auto_k)
fig, ax = plt.subplots(figsize=(6, 3.2))
rust_nmf.k_selection_plot(ax=ax)

# 2) Fit recommended recipe at auto_k
rust_nmf = ov.single.NMF(adata, rank=rust_nmf.auto_k, use_hvg=True, num_threads=8)
rust_nmf.fit(method='dnmf', init='nndsvd', max_iter=25)

# 3) Consensus heatmap (cNMF-style)
rust_nmf.consensus(n_runs=50, method='dnmf', max_iter=25)
rust_nmf.plot_consensus_heatmap(figsize=(6, 5))

# 4) Push results into adata
result_dict = rust_nmf.get_results(adata, key_added='NMF', n_top_genes=30)

# 5) RFC labels (cNMF-style)
rust_nmf.get_results_rfc(
    adata, result_dict,
    use_rep='scaled|original|X_pca',
    threshold=0.3, key_added='NMF_module_rfc',
)
```

## Constraints

- Do not pretend multi-worker `cNMF.factorize(...)` is complete unless every worker ran.
- Use `cNMF(..., use_gpu=True)` only when accelerator support is actually available and intended.
- Do not assume the notebook's chosen `k` or `density_threshold` transfers to other datasets.
- Do not run any RFC path without a real embedding in `obsm`.
- Treat `NMF.auto_k` as a recommendation, not a guarantee — always inspect `k_selection_plot` before committing.
- Do not use `ov.single.NMF` without `nmf-rs` installed; the wrapper raises lazily on first call.
- Keep smoke and acceptance commands shell-agnostic.

## Validation

**cNMF:**
- Check that `factorize(...)` actually ran the assigned iterations.
- Check that `consensus(...)` completed before `load_results(...)`.
- After `load_results(...)`, check that normalized usage columns exist.
- After `get_results(...)`, check that `cNMF_cluster` exists.
- After the RFC path, check that `cNMF_cluster_rfc` and `cNMF_cluster_clf` exist.

**NMF (Rust):**
- After `select_k_brunet(...)`, check `.auto_k` is populated and the returned DataFrame has one row per candidate K with silhouette + reconstruction columns.
- After `fit(...)`, check the factor matrices are populated (no NaNs).
- After `consensus(...)`, check `plot_consensus_heatmap` shows K block-diagonal structure; smeared edges → bump `n_runs` or switch `method='brunet'`.
- After `get_results(adata, key_added='NMF')`, check `adata.obsm['NMF_usage']`, `adata.varm['NMF_genes']`, `adata.obs['NMF_module']` all exist.
- After the RFC path, check `adata.obs[f'{key_added}_rfc']` exists.

- If only a bounded smoke path was run for either backend, say which expensive stages were reduced.

## Resource Map

- Use the branch selection notes when choosing CPU vs GPU, single-worker vs multi-worker, or direct labels vs RFC labels.
- Use the source grounding notes for current signatures, worker semantics, and output labels.
- Use the notebook mapping notes to trace the notebook's cNMF section into this reusable skill.
- Use the compatibility notes for compute and filesystem-sensitive behavior.
