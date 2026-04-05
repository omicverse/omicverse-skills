# Compatibility

- `Velo.dynamics(backend='scvelo')` uses scVelo recovery logic and can spawn worker processes; reviewer smoke should run from a real script context or force `fork` on macOS.
- `ov.pl.add_streamplot(...)` calls a grid helper that sets a local neighbor count from the number of cells. Tiny synthetic datasets can fail with `n_neighbors=0`; use at least 50 cells in smoke tests.
- `ov.pp.umap` and `ov.pp.leiden` depend on `ov.settings.mode`, so set that explicitly before validating backend behavior.
- `ov.pp.neighbors(method='rapids')` and the RAPIDS mode branches assume the RAPIDS stack is installed; keep the portable CPU path as the default smoke path.
- `Velo.cal_velocity(method='latentvelo')` is the heaviest branch and may require GPU or a representative import-only smoke instead of a full run.
- `Velo.cal_velocity(method='graphvelo')` is lighter than latentvelo but still branch-specific; keep reviewer data small and `n_jobs` low.
- `scv.datasets.dentategyrus()` is a tutorial convenience dataset and may download data; do not use it as the smoke default.
- `Velo.graphvelo(...)` can fail on tiny random gene subsets because its internal PCA expects enough features; for synthetic smoke, use a broad subset or all genes.
- `Velo.graphvelo(...)` writes projected vectors to `gv_X_umap` and `gv_X_pca` when `basis_keys=['X_umap', 'X_pca']`.
- The notebook's graph keys are `velocity_S_graph` and `velocity_gv_graph` for the respective layers; derive the embedding key from the chosen `velocity_key`.
