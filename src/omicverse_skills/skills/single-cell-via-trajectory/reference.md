# VIA trajectory — quick commands

## Vanilla VIA

```python
import omicverse as ov
import scanpy as sc
import matplotlib.pyplot as plt

ov.utils.ov_plot_set()

# Demo dataset (Setty 2019 hematopoiesis); replace with your AnnData
adata = ov.single.scRNA_hematopoiesis()
sc.tl.pca(adata, svd_solver='arpack', n_comps=200)

v0 = ov.single.pyVIA(
    adata=adata,
    adata_key='X_pca', adata_ncomps=80,
    basis='tsne',
    clusters='label',
    knn=30, random_seed=4,
    root_user=[4823],          # CD34+ HSC index in this cohort
)
v0.run()

# Topology
v0.plot_piechart_graph(clusters='label', cmap='Reds',
                        show_legend=False, ax_text=False, fontsize=4)

# Trajectory curves on the embedding
v0.plot_trajectory_gams(basis='tsne', clusters='label',
                         draw_all_curves=False)

# Streamplot
v0.plot_stream(basis='tsne', clusters='label',
               density_grid=0.8, scatter_size=30,
               scatter_alpha=0.3, linewidth=0.5)

# Streamplot coloured by pseudotime
v0.plot_stream(basis='tsne', density_grid=0.8, scatter_size=30,
               color_scheme='time', linewidth=0.5,
               min_mass=1, cutoff_perc=5,
               scatter_alpha=0.3, marker_edgewidth=0.1,
               density_stream=2, smooth_transition=1, smooth_grid=0.5)

# Lineage probabilities
v0.plot_lineage_probability(figsize=(8, 4))
v0.plot_lineage_probability(figsize=(6, 3), marker_lineages=[2, 3])

# Per-gene GAM trends
markers = ['IL3RA', 'IRF8', 'GATA1', 'GATA2', 'ITGA2B', 'MPO',
           'CD79B', 'SPI1', 'CD34', 'CSF1R', 'ITGAX']
v0.plot_gene_trend(gene_list=markers, figsize=(8, 6))
v0.plot_gene_trend_heatmap(gene_list=markers,
                            marker_lineages=[2], figsize=(4, 4))

# Cluster graph with gene overlays
v0.plot_clustergraph(gene_list=markers[:4], figsize=(12, 3))
```

## Velocity-guided VIA

```python
import omicverse as ov
import scanpy as sc
import scvelo as scv

ov.utils.ov_plot_set()

# 1) RNA velocity preprocessing (separate skill)
scv.pp.filter_and_normalize(adata, min_shared_counts=10, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)
sc.tl.pca(adata, n_comps=30)
sc.tl.umap(adata)

n_pcs = 30

# 2) Velocity-guided VIA
v0 = ov.single.pyVIA(
    adata=adata,
    adata_key='X_pca', adata_ncomps=n_pcs,
    basis='X_umap',
    clusters='clusters',
    knn=20, root_user=None,             # automatic root inference
    is_coarse=True,
    preserve_disconnected=True,
    pseudotime_threshold_TS=50,
    piegraph_arrow_head_width=0.15,
    piegraph_edgeweight_scalingfactor=2.5,
    velocity_matrix=adata.layers['velocity'],
    gene_matrix=adata.X.todense(),
    velo_weight=0.5,                    # 0=ignore velocity; 0.3-0.7 safe
    edgebundle_pruning_twice=False,
    edgebundle_pruning=0.15,
    pca_loadings=adata.varm['PCs'],
    random_seed=42,
)
v0.run()

# Same plotting helpers — basis is now UMAP-based
v0.plot_piechart_graph(clusters='clusters', cmap='Reds',
                        show_legend=False, ax_text=False, fontsize=4)
v0.plot_trajectory_gams(basis='X_umap', clusters='clusters',
                         draw_all_curves=False)
v0.plot_stream(basis='X_umap', clusters='clusters',
               density_grid=0.8, scatter_size=30,
               scatter_alpha=0.3, linewidth=0.5)
v0.plot_lineage_probability()
```

## Pseudotime extraction

```python
# Write pseudotime into adata.obs['via_pseudotime'] (or get as np.ndarray)
v0.get_pseudotime(adata)
print(adata.obs['via_pseudotime'].describe())
```

## Inspect lineage probabilities programmatically

```python
# Per-cluster lineage probabilities
piedict = v0.get_piechart_dict(label=0, clusters='label')
# piedict[cluster_id] -> {terminal_state: probability}
print(piedict)
```
