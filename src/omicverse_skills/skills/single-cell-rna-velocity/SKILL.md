---
name: omicverse-single-cell-rna-velocity
description: Analyze single-cell AnnData for RNA velocity with OmicVerse. Use when converting OmicVerse velocity notebooks into a reusable, triggerable skill, when deciding whether a velocity notebook subset or branch should update an existing skill, or when selecting the scvelo, dynamo, latentvelo, graphvelo, recipe, backend, or mode branches for velocity preprocessing, dynamics, and embedding.
---

# OmicVerse Single-Cell RNA Velocity

## Goal

Run reusable RNA velocity analysis on a velocity-ready `AnnData`: map spliced/unspliced layers if needed, fit moments and dynamics, estimate velocity, build the velocity graph, compute UMAP/Leiden, and render the stream plot. Use this same skill for method-comparison notebooks that branch across `scvelo`, `dynamo`, `latentvelo`, and `graphvelo`. If you still need FASTQ -> `h5ad` generation, hand off to the kb alignment skill first.

## Quick Workflow

1. Start from a velocity-ready `AnnData`.
2. Map source-specific layers such as `mature`/`nascent` to `spliced`/`unspliced` when needed.
3. Pick the branch trio up front: `recipe`, `backend`, and `method`.
4. Instantiate `ov.single.Velo(adata)` and run `filter_genes`, `preprocess`, `moments`, `dynamics`, and `cal_velocity`.
5. Build the velocity graph, neighbors, UMAP, Leiden, and velocity embedding.
6. If the notebook includes a GraphVelo refinement stage, run it after an initial velocity layer exists.
7. Plot the embedding and stream plot.
8. Validate the stored keys before treating the result as finished.

## Interface Summary

- `ov.single.Velo(adata)` creates the velocity workflow wrapper.
- `Velo.filter_genes(min_shared_counts=20)` keeps velocity-usable genes.
- `Velo.preprocess(recipe='monocle', n_neighbors=30, n_pcs=30, **kwargs)` uses Dynamo preprocessing and then builds a PCA-based neighbor graph.
- `Velo.moments(backend='dynamo'|'scvelo', ...)` and `Velo.dynamics(backend='dynamo'|'scvelo', ...)` branch on the requested backend.
- `Velo.cal_velocity(method='dynamo'|'scvelo'|'latentvelo'|'graphvelo', ...)` selects the velocity estimator and controls the output key family.
- `Velo.graphvelo(xkey='Ms', vkey='velocity_S', n_jobs=1, basis_keys=['X_umap', 'X_pca'], gene_subset=...)` refines an existing velocity layer and projects it to UMAP/PCA.
- `Velo.velocity_graph(vkey='velocity_S', ...)` and `Velo.velocity_embedding(basis='umap', vkey='velocity_S', ...)` project and store velocity outputs.
- `ov.pp.neighbors(..., method='umap'|'gauss'|'rapids', transformer=...)` builds the graph.
- `ov.pp.umap(adata, method='umap'|'rapids'|'torchdr'|'mde'|'pumap', **kwargs)` switches embedding backends.
- `ov.pp.leiden(...)` uses `ov.settings.mode` to choose CPU, mixed, or RAPIDS execution.
- `ov.pl.embedding(...)` and `ov.pl.add_streamplot(...)` render the final plot.

## Boundary

- Keep this skill on the analysis side of the notebook.
- Do not absorb FASTQ download, reference building, or kb counting here.
- If you need raw quantification, hand off to the kb alignment skill first.
- If the notebook only changes velocity `method`, `backend`, `recipe`, or `mode`, keep it in this skill instead of creating a duplicate skill.
- If the notebook introduces a truly different analysis family or input contract, split it.

## GraphVelo Refinement

- Use this branch after an initial velocity layer already exists in `adata`.
- Pick the matching input layer as `vkey` and use the corresponding `..._genes` mask when the notebook provides one.
- Call `Velo.graphvelo(..., basis_keys=['X_umap', 'X_pca'])` to project GraphVelo refinement back onto the existing manifold.
- Then call `Velo.velocity_graph(vkey='velocity_gv', xkey='Ms', n_jobs=8)` and `Velo.velocity_embedding(basis='umap', vkey='velocity_gv')`.
- Expect GraphVelo to add `velocity_gv`, `velocity_gv_genes`, `gv_X_umap`, and `gv_X_pca`; the downstream graph and embedding calls add `velocity_gv_graph` and `velocity_gv_umap`.

## Branch Selection

- Use `recipe='monocle'` for the notebook path; the live Dynamo `Preprocessor` also exposes `seurat`, `sctransform`, `pearson_residuals`, and `monocle_pearson_residuals`.
- Use `backend='scvelo'` for `moments` and `dynamics` when the input already has spliced/unspliced counts.
- Use `backend='dynamo'` when you want the Dynamo kinetics branch.
- Use `method='scvelo'` for `cal_velocity` on the portable notebook path.
- Use `method='dynamo'` when you want the Dynamo velocity branch.
- Use `method='latentvelo'` only when you have the required model inputs and a representative GPU or import-only smoke path.
- Use `method='graphvelo'` when you want the graph-based velocity branch; keep `n_jobs` small in smoke tests.
- Use `method='umap'` in `ov.pp.neighbors` unless you intentionally want the `gauss` or `rapids` graph branch.
- Set `ov.settings.mode` explicitly before `ov.pp.umap` and `ov.pp.leiden` if you need reproducible CPU, mixed, or GPU behavior.
- Use `ov.settings.cpu_init()` for the portable path; use `cpu_gpu_mixed_init()` or `gpu_init()` only when the hardware and dependencies actually support those modes.
- For notebook comparisons, keep the same input object and vary only the branch parameters you are studying.

## Input Contract

- Begin with an `AnnData` that already contains raw velocity counts or velocity-ready layers.
- If the source layers are named `mature` and `nascent`, rename or copy them to `spliced` and `unspliced`.
- `preprocess` should create `X_pca`; `neighbors` should then consume `X_pca`.
- `velocity_embedding` and the stream plot need `X_umap`.
- Keep the dataset large enough for stream plotting; tiny synthetic sets can fail because the grid helper derives its own neighbor count from cell number.

## Minimal Execution Patterns

```python
import omicverse as ov
import scanpy as sc

adata = sc.read_h5ad("velocity-ready.h5ad")
if "spliced" not in adata.layers and "mature" in adata.layers:
    adata.layers["spliced"] = adata.layers["mature"]
if "unspliced" not in adata.layers and "nascent" in adata.layers:
    adata.layers["unspliced"] = adata.layers["nascent"]

ov.settings.cpu_init()
velo = ov.single.Velo(adata)
velo.filter_genes(min_shared_counts=20)
velo.preprocess(recipe="monocle", n_neighbors=30, n_pcs=30)
velo.moments(backend="scvelo", n_neighbors=30, n_pcs=30)
velo.dynamics(backend="scvelo", n_jobs=1)
velo.cal_velocity(method="scvelo")
velo.velocity_graph(vkey="velocity_S", xkey="Ms", n_jobs=1)
ov.pp.neighbors(adata, n_neighbors=15, use_rep="X_pca")
ov.pp.umap(adata)
ov.pp.leiden(adata, resolution=0.2)
velo.velocity_embedding(basis="umap", vkey="velocity_S")
```

```python
# Method-comparison branch sketch
velo.moments(backend="dynamo")
velo.dynamics(backend="dynamo")
velo.cal_velocity(method="dynamo")

velo.cal_velocity(method="graphvelo", n_jobs=1)

velo.cal_velocity(
    method="latentvelo",
    batch_key="batch",
    celltype_key="celltype",
    velocity_key="velo_latentvelo",
)
```

## Validation

- Confirm `X_pca`, `X_umap`, and `leiden` are present after preprocessing and clustering.
- Confirm `Ms`, `Mu`, the chosen `velocity_key`, and the derived graph/embedding keys are present after the velocity steps.
- When you keep the default `velocity_key='velocity_S'`, confirm `velocity_S_graph` and `velocity_S_umap`.
- Confirm `ov.pl.add_streamplot(...)` returns an axes object and does not raise.
- If you switch `ov.settings.mode`, confirm the backend-specific messages match the selected mode.
- Validate `latentvelo` and any GPU-heavy path with a representative smoke or import-only interface check if a full run is not practical.
- For GraphVelo smoke on synthetic data, use a broad gene subset so the internal PCA has enough features.

## Resource Map

- Read `references/branch-selection.md` when choosing the recipe, backend, method, or mode branch.
- Read `references/source-notebook-map.md` when you need to see how the notebook splits into comparison branches and where the skill boundary sits.
- Read `references/source-grounding.md` for inspected signatures and source behavior.
- Read `references/compatibility.md` for smoke-path caveats and runtime notes.
