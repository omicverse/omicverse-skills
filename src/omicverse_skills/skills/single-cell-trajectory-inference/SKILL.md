---
name: omicverse-single-cell-trajectory-inference
description: Run or adapt OmicVerse single-cell trajectory inference on cluster-ready AnnData. Use when converting OmicVerse trajectory notebooks into a reusable skill, or when choosing the diffusion_map, slingshot, palantir, PAGA, or Palantir branch-selection branches for developmental ordering and lineage summaries.
---

# OmicVerse Single-Cell Trajectory Inference

## Goal

Run the reusable trajectory-analysis part of the OmicVerse notebook on a cluster-ready `AnnData`: choose a `TrajInfer` method branch, set origin and optional terminal states, compute pseudotime or fate outputs, and optionally summarize the trajectory with PAGA. Keep raw QC and preprocessing outside this skill; hand off to the preprocessing skill first when the object is not already graph-ready.

## Quick Workflow

1. Inspect the input `AnnData` and confirm it already has a usable representation, plotting basis, and cluster labels.
2. Choose the trajectory branch up front: `method='diffusion_map'`, `method='slingshot'`, or `method='palantir'`.
3. Instantiate `ov.single.TrajInfer(...)` with explicit `basis`, `use_rep`, `n_comps`, `n_neighbors`, and `groupby`.
4. Call `set_origin_cells(...)` before every trajectory run; call `set_terminal_cells(...)` when the branch supports or benefits from constrained terminal states.
5. Run `TrajInfer.inference(...)`.
6. If you need a coarse lineage graph, run `ov.utils.cal_paga(...)` on the resulting pseudotime and then `ov.utils.plot_paga(...)`.
7. If you chose `palantir`, optionally run `palantir_cal_branch(...)` and then `palantir_cal_gene_trends(...)` only when the extra dependencies and expression layer are available.
8. Validate the expected `obs`, `obsm`, `varm`, and `uns` keys before treating the result as finished.

## Interface Summary

- `ov.single.TrajInfer(adata, basis='X_umap', use_rep='X_pca', n_comps=50, n_neighbors=15, groupby='clusters')` creates the shared trajectory wrapper.
- `TrajInfer.set_origin_cells(origin)` stores the root label used by `diffusion_map`, `slingshot`, and `palantir`.
- `TrajInfer.set_terminal_cells(terminal)` stores optional terminal labels for `slingshot` and `palantir`.
- `TrajInfer.inference(method='palantir'|'diffusion_map'|'slingshot', **kwargs)` dispatches the main branch.
- `TrajInfer.palantir_cal_branch(q=..., eps=..., masks_key=...)` stores branch masks after a Palantir run.
- `TrajInfer.palantir_cal_gene_trends(layers='MAGIC_imputed_data')` computes lineage-wise trends after Palantir branch masks exist.
- `TrajInfer.palantir_plot_gene_trends(genes)` and `TrajInfer.palantir_plot_pseudotime(...)` are optional visualization helpers after the corresponding Palantir outputs exist.
- `ov.utils.cal_paga(adata, groups=..., vkey='paga', use_time_prior=..., minimum_spanning_tree=True)` computes the coarse topology summary.
- `ov.utils.plot_paga(...)` overlays the PAGA graph on the chosen basis or draws the graph layout directly.

Read `references/source-grounding.md` before adding more interface-specific details.

## Boundary

- Keep this skill on the trajectory-analysis side of the notebook.
- Do not absorb QC, HVG selection, PCA building, clustering, or general embedding generation here; those are already covered by the preprocessing skill.
- Keep `diffusion_map`, `slingshot`, `palantir`, PAGA, and Palantir follow-ups in one skill because they share the same cluster-ready `AnnData` contract and the same `TrajInfer` entrypoint.
- Split `sctour` into a separate skill because it changes the dependency and compute profile and is independently triggerable.

## Branch Selection

- Use `method='diffusion_map'` when you want DPT-style pseudotime from the same representation that drives the neighborhood graph.
- Use `method='slingshot'` when you need lineage curves on a low-dimensional embedding and can satisfy its extra dependency.
- Use `method='palantir'` when you need pseudotime plus entropy and fate probabilities, or when you plan to compute branch masks.
- Use `set_terminal_cells(...)` for Palantir when you know the terminal states and want to constrain the fate labels; leave them unset only when you intentionally want automatic terminal-state discovery.
- Run `ov.utils.cal_paga(...)` after you already have a pseudotime key such as `dpt_pseudotime` or `palantir_pseudotime`.
- Run `palantir_cal_branch(...)` only after `palantir_fate_probabilities` exists.
- Run `palantir_cal_gene_trends(...)` only after branch masks exist and you have a suitable expression layer.

## Input Contract

- Start from a cluster-ready `AnnData`.
- Ensure `groupby` points to a valid categorical key in `adata.obs`.
- Ensure `use_rep` points to a representation already stored in `adata.obsm`; the common portable path is `X_pca`.
- Ensure `basis` points to a plotting embedding already stored in `adata.obsm`; the notebook path uses `X_umap`.
- Ensure the origin label exists in `adata.obs[groupby]`.
- Ensure terminal labels exist in `adata.obs[groupby]` before calling `set_terminal_cells(...)`.
- Ensure a neighbor graph exists before `ov.utils.cal_paga(...)`.
- Ensure the expression layer used by `palantir_cal_gene_trends(...)` exists before calling it.

## Minimal Execution Patterns

```python
import omicverse as ov

traj = ov.single.TrajInfer(
    adata,
    basis="X_umap",
    use_rep="X_pca",
    n_comps=50,
    n_neighbors=15,
    groupby="clusters",
)
traj.set_origin_cells("Ductal")
traj.inference(method="diffusion_map")

ov.utils.cal_paga(
    adata,
    groups="clusters",
    vkey="paga",
    use_time_prior="dpt_pseudotime",
)
ov.utils.plot_paga(
    adata,
    basis="umap",
    color="clusters",
    title="PAGA DPT graph",
    min_edge_width=2,
    node_size_scale=1.5,
    show=False,
)
```

```python
import omicverse as ov

traj = ov.single.TrajInfer(
    adata,
    basis="X_umap",
    use_rep="X_pca",
    n_comps=50,
    n_neighbors=15,
    groupby="clusters",
)
traj.set_origin_cells("Ductal")
traj.set_terminal_cells(["Alpha", "Beta", "Delta", "Epsilon"])
traj.inference(method="palantir", num_waypoints=500)
traj.palantir_cal_branch(eps=0)

gene_trends = traj.palantir_cal_gene_trends(layers="MAGIC_imputed_data")
traj.palantir_plot_gene_trends(["Pax4", "Ins2"])
```

```python
import matplotlib.pyplot as plt
import omicverse as ov

fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(8, 8))
traj = ov.single.TrajInfer(
    adata,
    basis="X_umap",
    use_rep="X_pca",
    groupby="clusters",
)
traj.set_origin_cells("Ductal")
traj.inference(method="slingshot", num_epochs=1, debug_axes=axes)
```

## Constraints

- Do not hardcode notebook-specific cell types, terminal labels, or marker genes as if they were universal defaults.
- Do not rebuild preprocessing here when the real problem is that the input is not trajectory-ready; hand off to the preprocessing skill instead.
- Set `basis`, `use_rep`, `groupby`, and `method` explicitly; do not rely on notebook-local state.
- `slingshot` needs the extra `pcurvepy2` dependency in the current OmicVerse implementation.
- Palantir trend computation needs the extra `mellon` dependency.
- The current Palantir wrapper does not expose a direct `run_magic_imputation(..., n_jobs=...)` override; keep local smoke-specific workarounds in the validation harness, not in the reusable workflow.
- Keep notebook display-polish cells optional unless the user explicitly wants the same figure layout.

## Validation

- For `diffusion_map`, confirm `adata.obs['dpt_pseudotime']` exists after inference.
- For `slingshot`, confirm `adata.obs['slingshot_pseudotime']` exists after inference.
- For `palantir`, confirm `adata.obs['palantir_pseudotime']`, `adata.obs['palantir_entropy']`, and `adata.obsm['palantir_fate_probabilities']` exist after inference.
- After `palantir_cal_branch(...)`, confirm `adata.obsm['branch_masks']` exists.
- After `palantir_cal_gene_trends(...)`, confirm one or more `palantir_gene_trends_*` or `gene_trends_*` keys exist in `adata.varm`.
- After `ov.utils.cal_paga(...)`, confirm `adata.uns['paga']` exists with connectivity and transition-confidence entries.
- If you only validated a bounded smoke path, say so; do not claim full notebook reproduction.

## Resource Map

- Read `references/branch-selection.md` when deciding whether the notebook change belongs in this skill or the separate `sctour` skill, and when choosing a `TrajInfer` method branch.
- Read `references/source-notebook-map.md` to see which notebook cells were reused, delegated to the preprocessing skill, or split into the separate `sctour` skill.
- Read `references/source-grounding.md` for checked signatures, source-derived branch behavior, and dependency-sensitive follow-up steps.
- Read `references/compatibility.md` for dependency caveats and current runtime limitations.
- Run `scripts/trajectory_smoke.py` only for bounded local validation of representative branches; it is a smoke utility, not the reusable workflow contract.
