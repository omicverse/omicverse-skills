# Branch Selection

## Notebook Path

- `Velo.preprocess(recipe='monocle')`
- `Velo.moments(backend='scvelo')`
- `Velo.dynamics(backend='scvelo')`
- `Velo.cal_velocity(method='scvelo')`
- `ov.pp.neighbors(method='umap')`
- `ov.pp.umap(..., method='umap')`
- `ov.pp.leiden(...)` with `ov.settings.mode='cpu'`
- `Velo.graphvelo(xkey='Ms', vkey=<existing velocity layer>, basis_keys=['X_umap', 'X_pca'])`
- `Velo.velocity_graph(vkey='velocity_gv', xkey='Ms')`
- `Velo.velocity_embedding(basis='umap', vkey='velocity_gv')`

## Live Branches Worth Remembering

- `recipe` branches from Dynamo preprocessing: `monocle`, `seurat`, `sctransform`, `pearson_residuals`, `monocle_pearson_residuals`.
- `backend` branches in velocity moments/dynamics: `dynamo` and `scvelo`.
- `method` branches in velocity estimation: `dynamo`, `scvelo`, `latentvelo`, `graphvelo`.
- `method` branches in neighbor search: `umap`, `gauss`, `rapids`.
- `method` branches in UMAP embedding: `umap`, `rapids`, `torchdr`, `mde`, `pumap`.
- `ov.settings.mode` branches in preprocessing UMAP/Leiden: `cpu`, `cpu-gpu-mixed`, and GPU/RAPIDS mode.
- `graphvelo` takes an existing velocity layer plus `basis_keys`; the notebook uses `['X_umap', 'X_pca']`.

## Selection Rules

- Use `scvelo` when the input already has spliced/unspliced counts and you want the portable notebook path.
- Use `dynamo` when you need the Dynamo kinetics branch or the notebook is comparing estimators.
- Use `latentvelo` only when you have the required cell metadata and a representative GPU or import-only smoke path.
- Use `graphvelo` when you need the downstream GraphVelo refinement branch and already have a velocity layer to refine.
- Use `graphvelo` on a broad gene subset in synthetic smoke; tiny random subsets can underflow the internal PCA rank.
- Use `cpu` for the portable smoke path.
- Use `cpu-gpu-mixed` only when you actually want the torch-backed mixed path.
- Use `rapids` only if the RAPIDS stack is installed and you want the GPU backend.
- For notebook comparisons, keep the shared preprocessing stable and vary only the studied branch parameter.
