# Source Grounding Notes

## Core integration entrypoint

`ov.single.batch_correction(...)`

Observed signature:

```python
(adata: anndata._core.anndata.AnnData, batch_key: str, use_rep='scaled|original|X_pca', methods: str = 'harmony', n_pcs: int = 50, **kwargs)
```

Grounded method branches in live source:

- `harmony`
- `combat`
- `scanorama`
- `scVI`
- `CellANOVA`
- `Concord`
- `concord`

## Backend-specific behavior

### Harmony

- If `use_rep='scaled|original|X_pca'` is requested and that key is missing, the function scales and runs PCA internally before integration.
- Live source calls `run_harmony(...)` and writes:
  - `adata.obsm['X_pca_harmony']`
  - `adata.obsm['X_harmony']`

`run_harmony(...)` signature:

```python
(data_mat: numpy.ndarray, meta_data: pandas.core.frame.DataFrame, vars_use, theta=None, lamb=None, sigma=0.1, nclust=None, tau=0, block_size=0.05, max_iter_harmony=10, max_iter_kmeans=20, epsilon_cluster=1e-05, epsilon_harmony=0.0001, plot_convergence=False, verbose=True, reference_values=None, cluster_prior=None, random_state=0, cluster_fn='kmeans', use_gpu=True, **kwargs)
```

Important branch knobs exposed here include `use_gpu`, `cluster_fn`, `theta`, `lamb`, and convergence controls.

### Combat

- Live source runs `scanpy.pp.combat(...)` on a copy, then rescales and reruns PCA.
- It writes `adata.obsm['X_combat']`.
- This branch returns `adata`, not a model object.

### Scanorama

- Live source expects the batch column to have categorical categories and iterates through those categories to split the input.
- It calls `integrate_scanpy(adatas, **kwargs)` and then concatenates the per-batch `X_scanorama` blocks.
- It writes `adata.obsm['X_scanorama']`.

`integrate_scanpy(...)` signature:

```python
(adatas, **kwargs)
```

### scVI

- Live source calls `scvi.model.SCVI.setup_anndata(adata, layer='counts', batch_key=batch_key)` before model construction.
- The branch writes `adata.obsm['X_scVI']`.
- The branch returns the trained scVI model object, not `adata`.

`scvi.model.SCVI(...)` signature:

```python
(adata: 'AnnData | None' = None, n_hidden: 'int' = 128, n_latent: 'int' = 10, n_layers: 'int' = 1, dropout_rate: 'float' = 0.1, dispersion: "Literal['gene', 'gene-batch', 'gene-label', 'gene-cell']" = 'gene', gene_likelihood: "Literal['zinb', 'nb', 'poisson', 'normal']" = 'zinb', latent_distribution: "Literal['normal', 'ln']" = 'normal', **kwargs)
```

Notebook-specific kwargs such as `n_layers`, `n_latent`, and `gene_likelihood='nb'` match the live constructor.

### CellANOVA

- Live source normalizes the HVG column naming mismatch by copying `highly_variable_features` into `highly_variable` when needed.
- It runs `calc_ME(...)`, then `calc_BE(...)`, then `calc_TE(...)`.
- It converts `adata.layers['denoised']` to CSR, runs PCA on `layer='denoised'`, and writes `adata.obsm['X_cellanova']` from `denoised|original|X_pca`.

`calc_BE(...)` signature:

```python
(adata, integrate_key, control_dict, var_cutoff=0.9, k_max=1500, verbose=False, k_select=None)
```

`control_dict` is therefore not notebook sugar; it is a true required branch input for this path.

### Concord

- Live source accepts both `methods='Concord'` and `methods='concord'`.
- It chooses features from `highly_variable` first, then `highly_variable_features`, then falls back to all genes.
- It constructs `concord.Concord(..., domain_key=batch_key, **kwargs)`.
- It then calls `fit_transform(output_key='X_concord')`.
- The branch returns the Concord object.

`concord.Concord(...)` signature:

```python
(adata, save_dir='save/', copy_adata=False, verbose=False, **kwargs)
```

`concord.Concord.fit_transform(...)` signature:

```python
(self, output_key='Concord', return_decoded=False, decoder_domain=None, return_class=True, return_class_prob=True, save_model=True)
```

## Visualization helpers

Observed signatures:

```python
ov.pp.mde(adata, embedding_dim=2, n_neighbors=15, basis='X_mde', n_pcs=None, use_rep=None, knn=True, transformer=None, metric='euclidean', verbose=False, key_added=None, random_state=0, repulsive_fraction=0.7, constraint=None)
```

```python
ov.utils.mde(data, device=None, **kwargs)
```

The notebook uses `ov.utils.mde(...)` on integrated embeddings and stores the result manually into `adata.obsm[...]`.

## Benchmarking entrypoint

`scib_metrics.benchmark.Benchmarker(...)`

Observed signature:

```python
(adata: anndata._core.anndata.AnnData, batch_key: str, label_key: str, embedding_obsm_keys: list[str], bio_conservation_metrics: scib_metrics.benchmark._core.BioConservation | None = BioConservation(isolated_labels=True, nmi_ari_cluster_labels_leiden=False, nmi_ari_cluster_labels_kmeans=True, silhouette_label=True, clisi_knn=True), batch_correction_metrics: scib_metrics.benchmark._core.BatchCorrection | None = BatchCorrection(silhouette_batch=True, ilisi_knn=True, kbet_per_label=True, graph_connectivity=True, pcr_comparison=True), pre_integrated_embedding_obsm_key: str | None = None, n_jobs: int = 1, progress_bar: bool = True)
```

Method signatures:

```python
prepare(self, neighbor_computer=None)
benchmark(self)
plot_results_table(self, min_max_scale=True, show=True, save_dir=None)
```

Grounded benchmark behavior from live source:

- If `pre_integrated_embedding_obsm_key` is `None`, `prepare()` computes `X_pca` itself with `scanpy.tl.pca(..., use_highly_variable=False)`.
- Each embedding key becomes a temporary `AnnData` with copied batch and label columns.
- Neighbor graphs are prepared at `15`, `50`, and `90` neighbors.
- The benchmark results table is indexed by metric name and includes one column per embedding key.
