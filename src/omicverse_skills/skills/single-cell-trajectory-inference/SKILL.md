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
6. For lineage topology, use the unified `ov.pl.trajectory(...)` / `ov.pl.trajectory_overlay(...)` plotters — they consume the trajectory state written into `adata` and share the same visual grammar across `diffusion_map` / `slingshot` / `palantir` / `monocle` / `sctour` / CellRank. `ov.utils.cal_paga(...)` + `ov.utils.plot_paga(...)` remain valid for explicit PAGA computation.
7. Summarize the branch structure with `ov.pl.branch_streamplot(adata, group_key=..., pseudotime_key=...)` — river-style plot driven only by a pseudotime vector and cluster labels, interoperable across methods.
8. Fit marker trends with the shared GAM stack: `ov.single.dynamic_features(adata, genes=..., pseudotime=...)` → `ov.pl.dynamic_trends(res, ...)` for per-gene curves and `ov.pl.dynamic_heatmap(adata, var_names=..., pseudotime=...)` for many-gene panels. The same backend is used by every trajectory skill, so output is interoperable.
9. If you chose `palantir`, optionally run `palantir_cal_branch(...)` and then `palantir_cal_gene_trends(...)` only when the extra dependencies and expression layer are available.
10. Validate the expected `obs`, `obsm`, `varm`, and `uns` keys before treating the result as finished.

## Interface Summary

**Trajectory fitting (`TrajInfer`)**

- `ov.single.TrajInfer(adata, basis='X_umap', use_rep='X_pca', n_comps=50, n_neighbors=15, groupby='clusters')` creates the shared trajectory wrapper.
- `TrajInfer.set_origin_cells(origin)` stores the root label used by `diffusion_map`, `slingshot`, and `palantir`.
- `TrajInfer.set_terminal_cells(terminal)` stores optional terminal labels for `slingshot` and `palantir`.
- `TrajInfer.inference(method='palantir'|'diffusion_map'|'slingshot', **kwargs)` dispatches the main branch.
- `TrajInfer.palantir_cal_branch(q=..., eps=..., masks_key=...)` stores branch masks after a Palantir run.
- `TrajInfer.palantir_cal_gene_trends(layers='MAGIC_imputed_data')` computes lineage-wise trends after Palantir branch masks exist.
- `TrajInfer.palantir_plot_gene_trends(genes)` and `TrajInfer.palantir_plot_pseudotime(...)` are optional visualization helpers after the corresponding Palantir outputs exist.

**Unified `ov.pl` trajectory plotting (shared across methods — added by commit `4f28ab6`)**

- `ov.pl.trajectory(adata, *, method='diffusion'|'slingshot'|'palantir'|'monocle'|'sctour'|'cellrank', basis='X_umap', color=...)` draws the trajectory backbone on the chosen embedding. Replaces method-specific plot calls.
- `ov.pl.trajectory_overlay(adata, *, ax, method=..., ...)` overlays the trajectory backbone on an existing axis already populated by `ov.pl.embedding(...)`. Use when you want a custom-styled UMAP with the PAGA / curve overlay on top.
- `ov.pl.branch_streamplot(adata, *, group_key, pseudotime_key, trunk_groups=None, branch_center=0.5, figsize=..., show=False)` — river-style per-cluster pseudotime view. Works with any pseudotime key (`dpt_pseudotime`, `palantir_pseudotime`, `slingshot_pseudotime`, `sctour_pseudotime`, `velocity_pseudotime`, etc.).
- `ov.pl.dynamic_heatmap(adata, *, pseudotime, var_names, cell_annotation=None, use_fitted=True, cell_bins=200, ...)` — GAM-fitted heatmap of many genes along pseudotime; `var_names` accepts a list or `{program_name: [genes]}` dict for multi-module panels.
- `ov.utils.cal_paga(adata, groups=..., vkey='paga', use_time_prior=..., minimum_spanning_tree=True)` still computes an explicit PAGA topology summary; `ov.utils.plot_paga(...)` is still valid but the `ov.pl.trajectory*` family is the canonical entrypoint in the current tutorials.

**Marker-dynamics GAM stack (shared backend across all trajectory skills)**

- `ov.single.dynamic_features(adata, genes, pseudotime, *, layer=None, groupby=None, groups=None, n_splines=8, store_raw=True, raw_obs_keys=..., key_added='dynamic_features')` fits per-gene GAM curves. Pass `groupby + groups` for branch-aware fits.
- `ov.pl.dynamic_trends(res, *, genes, compare_features=False, compare_groups=False, split_time=None, shared_trunk=True, add_point=True, point_color_by=..., line_style_by=..., ...)` plots single-line global trends, multi-marker overlays, or branch-aware comparisons.

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

# Unified trajectory plotting — canonical in the current tutorials
ov.pl.trajectory(adata, method="diffusion", basis="X_umap", color="clusters")

# Or overlay on an existing custom embedding
fig, ax = ov.plt.subplots(figsize=(4, 4))
ov.pl.embedding(adata, basis="X_umap", color="clusters", ax=ax, show=False)
ov.pl.trajectory_overlay(adata, ax=ax, method="diffusion")

# River-style branch view
ov.pl.branch_streamplot(
    adata, group_key="clusters", pseudotime_key="dpt_pseudotime", show=False,
)

# Marker trends along DPT pseudotime
res = ov.single.dynamic_features(
    adata, genes=["Pdx1", "Ins1", "Gcg"],
    pseudotime="dpt_pseudotime", store_raw=True, raw_obs_keys=["clusters"],
)
ov.pl.dynamic_trends(
    res, genes=["Pdx1", "Ins1", "Gcg"],
    add_point=True, point_color_by="clusters",
)

# Optional: explicit PAGA still supported
ov.utils.cal_paga(adata, groups="clusters", vkey="paga", use_time_prior="dpt_pseudotime")
ov.utils.plot_paga(adata, basis="umap", color="clusters", show=False)
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
- After `ov.pl.trajectory(...)` / `ov.pl.trajectory_overlay(...)`, the figure should show the backbone connecting clusters in a topology consistent with the chosen `method`; an empty overlay means the trajectory state was not written into `adata` (re-run `TrajInfer.inference(...)` or `ov.utils.cal_paga(...)`).
- After `ov.single.dynamic_features(...)`, the returned result should have `.failed` (list of skipped genes with reasons) and a populated per-gene fit table; treat `.failed` non-empty as expected for low-variance / sparse genes.
- If you only validated a bounded smoke path, say so; do not claim full notebook reproduction.

## Resource Map

- Read `references/branch-selection.md` when deciding whether the notebook change belongs in this skill or the separate `sctour` skill, and when choosing a `TrajInfer` method branch.
- Read `references/source-notebook-map.md` to see which notebook cells were reused, delegated to the preprocessing skill, or split into the separate `sctour` skill.
- Read `references/source-grounding.md` for checked signatures, source-derived branch behavior, and dependency-sensitive follow-up steps.
- Read `references/compatibility.md` for dependency caveats and current runtime limitations.
- Run `scripts/trajectory_smoke.py` only for bounded local validation of representative branches; it is a smoke utility, not the reusable workflow contract.
