---
name: omicverse-single-cell-sctour-trajectory
description: Run the OmicVerse sctour trajectory branch on raw-count single-cell AnnData. Use when adapting the scTour part of an OmicVerse trajectory notebook, or when you need sctour pseudotime, latent space, or vector-field outputs instead of the diffusion_map, slingshot, or palantir branches.
---

# OmicVerse Single-Cell scTour Trajectory

## Goal

Run the `sctour` branch of `ov.single.TrajInfer` on a raw-count `AnnData` and capture the trainer-produced pseudotime, latent space, and vector field. Keep this separate from the lighter graph-based trajectory skill because it has a different dependency and compute profile.

## Quick Workflow

1. Inspect the input `AnnData` and confirm raw UMI counts are available in `.X`.
2. Preserve or restore counts into `.X` before the `sctour` call.
3. Instantiate `ov.single.TrajInfer(...)` with the plotting `basis`, `groupby`, and any existing representation you want to preserve for downstream plots.
4. Run `TrajInfer.inference(method='sctour', **trainer_kwargs)`.
5. Validate that pseudotime, latent coordinates, and vector-field outputs were written.
6. For lineage topology, draw the unified trajectory backbone with `ov.pl.trajectory(adata, method='sctour', ...)` or overlay on a custom UMAP with `ov.pl.trajectory_overlay(adata, ax=ax, method='sctour')` — same visual grammar shared with the other trajectory skills (commit `4f28ab6`).
7. Summarize branch structure with `ov.pl.branch_streamplot(adata, group_key=..., pseudotime_key='sctour_pseudotime')`.
8. Fit marker trends with the shared GAM stack: `ov.single.dynamic_features(adata, genes=..., pseudotime='sctour_pseudotime')` → `ov.pl.dynamic_trends(...)`; many-gene panels via `ov.pl.dynamic_heatmap(...)`.
9. Only invert or rescale pseudotime after fit when that direction change is biologically justified for the dataset.

## Interface Summary

**Fitting**

- `ov.single.TrajInfer(...)` is still the public entrypoint.
- `TrajInfer.inference(method='sctour', **kwargs)` dispatches the `sctour` branch.
- In the current wrapper source, the branch internally calls `sct.train.Trainer(self.adata, loss_mode='nb', **kwargs)`.
- The wrapper then calls `tnode.train()`, stores `adata.obs['sctour_pseudotime'] = tnode.get_time()`, stores `adata.obsm['X_TNODE']` from `tnode.get_latentsp(alpha_z=0.5, alpha_predz=0.5)`, and stores `adata.obsm['X_VF']` from `tnode.get_vector_field(...)`.

**Visualisation (shared with other trajectory skills)**

- `ov.pl.trajectory(adata, method='sctour', basis='X_umap', color=...)` — unified trajectory plot.
- `ov.pl.trajectory_overlay(adata, ax=ax, method='sctour')` — overlay backbone on existing embedding.
- `ov.pl.branch_streamplot(adata, group_key=..., pseudotime_key='sctour_pseudotime', show=False)` — river-style view.
- `ov.single.dynamic_features(adata, genes=..., pseudotime='sctour_pseudotime', ...)` + `ov.pl.dynamic_trends(...)` — per-gene GAM curves.
- `ov.pl.dynamic_heatmap(adata, pseudotime='sctour_pseudotime', var_names=..., ...)` — many-gene panel.

Read `references/source-grounding.md` before adding more interface-specific detail.

## Boundary

- Keep `sctour` separate from the shared `diffusion_map` / `slingshot` / `palantir` skill.
- Do not fold preprocessing or clustering into this skill unless the real request is upstream data preparation; hand off to the preprocessing skill first.
- Do not use this skill when the user only wants DPT pseudotime, lineage curves, or Palantir fate probabilities; use the trajectory-inference skill instead.

## Stage Selection

- Use this skill only when the requested branch is explicitly `method='sctour'` or when the user asks for scTour latent-space or vector-field outputs.
- Keep the notebook's `alpha_recon_lec` and `alpha_recon_lode` values as branch examples, not universal defaults.
- Treat the wrapper's fixed `loss_mode='nb'` as part of the current execution path; if a different loss mode is needed, inspect the underlying `sctour` package before changing the workflow.

## Input Contract

- Start from `AnnData`.
- Ensure raw UMI counts are present in `.X` for the default negative-binomial loss path.
- Preserve a plotting basis such as `X_umap` only if you want to overlay `sctour_pseudotime` later.
- Preserve cluster annotations only if you need grouped plots or interpretation after training.

## Minimal Execution Pattern

```python
import omicverse as ov

traj = ov.single.TrajInfer(
    adata,
    basis="X_umap",
    use_rep="X_pca",
    groupby="clusters",
)
traj.inference(
    method="sctour",
    alpha_recon_lec=0.5,
    alpha_recon_lode=0.5,
)

# Optional notebook-specific direction flip only if biologically justified.
adata.obs["sctour_pseudotime"] = 1 - adata.obs["sctour_pseudotime"]

# Unified trajectory plots
ov.pl.trajectory(adata, method="sctour", basis="X_umap", color="clusters")

fig, ax = ov.plt.subplots(figsize=(4, 4))
ov.pl.embedding(adata, basis="X_umap", color="clusters", ax=ax, show=False)
ov.pl.trajectory_overlay(adata, ax=ax, method="sctour")

# Branch streamplot + marker dynamics on sctour_pseudotime
ov.pl.branch_streamplot(
    adata, group_key="clusters", pseudotime_key="sctour_pseudotime", show=False,
)
res = ov.single.dynamic_features(
    adata, genes=["Pdx1", "Ins1", "Gcg"],
    pseudotime="sctour_pseudotime", store_raw=True, raw_obs_keys=["clusters"],
)
ov.pl.dynamic_trends(res, genes=["Pdx1", "Ins1", "Gcg"],
                    add_point=True, point_color_by="clusters")
```

## Constraints

- Keep raw counts in `.X` before training; the wrapper hardcodes `loss_mode='nb'`.
- Treat the notebook's pseudotime inversion as optional and dataset-specific.
- This branch depends on the external `sctour` package and was not fully executable in the current local environment.
- Do not promise GPU use or full-scale training unless the runtime and dependency stack are confirmed.

## Validation

- Confirm `adata.obs['sctour_pseudotime']` exists after the fit.
- Confirm `adata.obsm['X_TNODE']` exists after the fit.
- Confirm `adata.obsm['X_VF']` exists after the fit.
- If you invert pseudotime, record that postprocessing step explicitly.
- If the environment cannot import `sctour`, treat the branch as source-grounded but not empirically reproduced.

## Resource Map

- Read `references/source-notebook-map.md` to see why `sctour` was split out from the shared trajectory skill.
- Read `references/source-grounding.md` for the current wrapper behavior and stored output keys.
- Read `references/compatibility.md` for dependency and validation limitations.
