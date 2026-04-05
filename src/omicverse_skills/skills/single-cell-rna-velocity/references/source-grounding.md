# Source Grounding

## Inspected Interfaces

- `omicverse.single.Velo(adata)`
- `Velo.filter_genes(min_shared_counts=20)`
- `Velo.preprocess(recipe='monocle', n_neighbors=30, n_pcs=30, **kwargs)`
- `Velo.moments(backend='dynamo', n_pcs=30, n_neighbors=30, **kwargs)`
- `Velo.dynamics(backend='dynamo', **kwargs)`
- `Velo.cal_velocity(method='dynamo', batch_key=None, celltype_key=None, velocity_key='velocity_S', n_jobs=1, n_top_genes=2000, param_name_key='tmp/latentvelo_params', latentvelo_VAE_kwargs={}, **kwargs)`
- `Velo.graphvelo(xkey='Ms', vkey='velocity_S', n_jobs=1, basis_keys=['X_umap', 'X_pca'], gene_subset=None, **kwargs)`
- `Velo.velocity_graph(basis='umap', vkey='velocity_S', **kwargs)`
- `Velo.velocity_embedding(basis='umap', vkey='velocity_S', **kwargs)`
- `ov.pp.neighbors(..., method='umap'|'gauss'|'rapids', transformer=...)`
- `ov.pp.umap(adata, method='umap'|'rapids'|'torchdr'|'mde'|'pumap', **kwargs)`
- `ov.pp.leiden(adata, resolution=1.0, random_state=0, key_added='leiden', local_iterations=100, max_levels=10, device='cpu', symmetrize=None, **kwargs)`
- `dynamo.preprocessing.Preprocessor.preprocess_adata(recipe='monocle'|'seurat'|'sctransform'|'pearson_residuals'|'monocle_pearson_residuals')`

## Source Notes

- `omicverse/single/_velo.py` defines the velocity wrapper methods and the `method`/`backend` branches used by the notebook.
- `omicverse/single/_velo.py` also defines `graphvelo`, which refines an existing velocity layer and writes `velocity_gv` plus projected basis-specific vectors.
- `omicverse/pp/_preprocess.py` routes `leiden` on `ov.settings.mode` and uses CPU, mixed, or RAPIDS paths.
- `omicverse/pp/_umap.py` exposes additional UMAP methods: `umap`, `torchdr`, `mde`, `pumap`, and `rapids`.
- `omicverse/pp/_neighbors.py` exposes `method='umap'|'gauss'|'rapids'` plus transformer branches such as `pynndescent`, `rapids`, `sklearn`, and `pyg`.
- `dynamo/preprocessing/Preprocessor.py` is the source of the recipe branch choices used by `Velo.preprocess`.

## Notebook-Grounded Claims

- The notebook's first branch uses `recipe='monocle'`, `backend='scvelo'`, and `method='scvelo'`.
- The notebook's second branch uses `method='dynamo'`.
- The notebook's third branch uses `method='latentvelo'` and a GPU-oriented configuration.
- The notebook's fourth branch uses `Velo.graphvelo(...)` followed by `velocity_graph(vkey='velocity_gv')` and `velocity_embedding(vkey='velocity_gv')`.
- The notebook's graph/embedding/plot path uses `ov.pp.neighbors`, `ov.pp.umap`, `ov.pp.leiden`, `Velo.velocity_embedding`, `ov.pl.embedding`, and `ov.pl.add_streamplot`.
