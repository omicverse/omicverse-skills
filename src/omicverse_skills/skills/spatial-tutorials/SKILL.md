---
name: spatial-transcriptomics-tutorials-with-omicverse
title: Spatial transcriptomics tutorials with omicverse
description: "Spatial transcriptomics: Visium/HD, Stereo-seq, Slide-seq preprocessing (crop, rotate, cellpose), deconvolution (Tangram, cell2location, Starfysh), clustering (GraphST, STAGATE), integration, trajectory, communication."
---

# Spatial Transcriptomics with OmicVerse

This skill covers spatial analysis workflows organized into three stages: Preprocessing, Deconvolution, and Downstream Analysis. Each stage includes the critical function calls, parameter guidance, and common pitfalls.

## Defensive Validation: Always Check Spatial Coordinates First

Before ANY spatial operation, verify that spatial coordinates exist and are numeric:

```python
# Required check before spatial analysis
assert 'spatial' in adata.obsm, \
    "Missing adata.obsm['spatial']. Load with ov.io.spatial.read_visium() or set manually."
# Cast to float64 to prevent coordinate precision issues during rotation/cropping
adata.obsm['spatial'] = adata.obsm['spatial'].astype('float64')
```

## Stage 1: Preprocessing

### Crop, Rotate, and Align Coordinates

Load Visium data and manipulate spatial coordinates for region selection and alignment:

```python
import scanpy as sc, omicverse as ov
ov.plot_set()

adata = sc.datasets.visium_sge(sample_id="V1_Breast_Cancer_Block_A_Section_1")
library_id = list(adata.uns['spatial'].keys())[0]

# Cast coordinates before manipulation
adata.obsm['spatial'] = adata.obsm['spatial'].astype('float64')

# Crop to region of interest
adata_crop = ov.space.crop_space_visium(adata, crop_loc=(0, 0), crop_area=(1000, 1000),
                                         library_id=library_id, scale=1)

# Rotate and auto-align
adata_rot = ov.space.rotate_space_visium(adata, angle=45, library_id=library_id)
ov.space.map_spatial_auto(adata_rot, method='phase')
# For manual refinement: ov.space.map_spatial_manual(adata_rot, ...)
```

### Visium HD Cell Segmentation

Segment Visium HD bins into cells using cellpose:

```python
adata = ov.space.read_visium_10x(path="binned_outputs/square_002um/",
                                  source_image_path="tissue_image.btf")
ov.pp.filter_genes(adata, min_cells=3)
ov.pp.filter_cells(adata, min_counts=1)

# H&E-based segmentation
adata = ov.space.visium_10x_hd_cellpose_he(adata, mpp=0.3, gpu=True, buffer=150)
# Expand labels to neighboring bins
ov.space.visium_10x_hd_cellpose_expand(adata, labels_key='labels_he',
                                        expanded_labels_key='labels_he_expanded', max_bin_distance=4)
# Gene-expression-driven seeds
ov.space.visium_10x_hd_cellpose_gex(adata, obs_key="n_counts_adjusted", mpp=0.3, sigma=5)
# Merge labels and aggregate to cell-level
ov.space.salvage_secondary_labels(adata, primary_label='labels_he_expanded',
                                   secondary_label='labels_gex', labels_key='labels_joint')
cdata = ov.space.bin2cell(adata, labels_key='labels_joint')
```

### Xenium Preprocessing

10x Genomics Xenium output (cell × gene matrices with polygon segmentation already done by the instrument). Read with `ov.io.read_xenium`; cache the parsed AnnData on disk for fast reuse:

```python
import omicverse as ov
ov.style(font_path='Arial')
ov.settings.cpu_gpu_mixed_init()        # optional: enable mixed CPU/GPU acceleration

# Load — set load_image=False to skip morphology image when not needed
adata = ov.io.read_xenium('data/xenium_breast_rep1', load_image=False)

# Cache for fast re-load (cold parse → write cache; warm reads are 10–50x faster)
adata = ov.io.read_xenium(
    'data/xenium_breast_rep1',
    cache_file='data/xenium_breast_rep1_cache.h5ad',
)

# Standard scanpy preprocessing
ov.pp.normalize_total(adata, target_sum=1e4)
ov.pp.log1p(adata)
ov.pp.scale(adata)
ov.pp.pca(adata, layer='scaled', n_pcs=50)
ov.pp.neighbors(adata, n_neighbors=15, use_rep='scaled|original|X_pca', n_pcs=50)
ov.pp.leiden(adata, resolution=0.5)

# Spatial visualisation — Xenium coordinates have inverted y; remember to invert axis
ov.pl.embedding(adata, basis='spatial', color='leiden',
                palette=ov.pl.palette_112, legend_fontsize=8)
```

For cell-segmentation overlay (polygon-aware figures with optional H&E / DAPI background), use `ov.pl.spatialseg`:

```python
library_id = list(adata.uns['spatial'].keys())[0]

ov.pl.spatialseg(
    adata, color='leiden',
    library_id=library_id,
    edges_color='white', edges_width=0.3,
    alpha=1.0, legend_fontsize=8,
    palette=ov.pl.palette_112,
    crop_coord=(2000, 3200, 2500, 3700),    # x0, x1, y0, y1 in spatial coords
    figsize=(7, 6),
)

# With morphology image overlay (DAPI / H&E behind cells)
adata_img = ov.io.read_xenium('data/xenium_breast_rep1', load_image=True)
# Re-attach the same processed obs / obsm before plotting
ov.pl.spatialseg(
    adata_img, color='KRT7',
    library_id=library_id,
    edges_color='white', edges_width=0.4,
    alpha=0.65, alpha_img=1.0,        # 0.45–0.65 keeps morphology visible
    cmap=ov.pl.create_custom_colormap('#a51616'), vmax=10,
    seg_contourpx=1.5,                # dashed cell-outline contour
    crop_coord=(2000, 3200, 2500, 3700),
    figsize=(7, 6),
)
```

Xenium-specific gotchas:
- **`vmax='p99.2'`** on `ov.pl.embedding` clips the long-tailed expression distribution that's typical of Xenium probes; without it a few hot cells dominate the colour scale.
- **`alpha=0.45–0.65`** on `spatialseg` plots when `alpha_img=1.0` keeps DAPI morphology readable through cluster fills. Higher alpha hides the background.
- **Spatial axis**: Xenium spatial coordinates use image-pixel convention (y increases downward); call `ax.invert_yaxis()` on `ov.pl.embedding` to match the morphology image orientation.
- Cache the parsed AnnData with `cache_file=` — the cold parse of a 100k-cell Xenium run is 10–60 s; warm reads are sub-second.

## Stage 2: Deconvolution

### Critical API Reference: Method Selection

```python
# CORRECT — Tangram: method passed to deconvolution() call
decov_obj = ov.space.Deconvolution(adata_sc=sc_adata, adata_sp=sp_adata,
                                    celltype_key='Subset', result_dir='result/tangram')
decov_obj.preprocess_sc(max_cells=5000)
decov_obj.preprocess_sp()
decov_obj.deconvolution(method='tangram', num_epochs=1000)

# CORRECT — cell2location: method passed at INIT time
cell2_obj = ov.space.Deconvolution(adata_sc=sc_adata, adata_sp=sp_adata,
                                    celltype_key='Subset', result_dir='result/c2l',
                                    method='cell2location')
cell2_obj.deconvolution(max_epochs=30000)
cell2_obj.save_model('result/c2l/model')
```

Note: For cell2location, the `method` parameter is set at initialization, not at the `deconvolution()` call. For Tangram, it's passed to `deconvolution()`.

### FlashDeconv — Atlas-Scale Sketching-Based Deconvolution

For Visium HD / Slide-seq / Stereo-seq (10⁵–10⁶ spots) where Tangram and cell2location become impractical, use `method='FlashDeconv'`. It runs on CPU, finishes in minutes for 10k spots vs. 60+ for cell2location, and ships built-in spatial regularisation:

```python
decov_obj = ov.space.Deconvolution(
    adata_sc=sc_adata,
    adata_sp=sp_adata,
    celltype_key='Subset',
    result_dir='result/flashdeconv',
    method='FlashDeconv',                 # set at INIT time, like cell2location
)
decov_obj.preprocess_sc(max_cells=5000)
decov_obj.preprocess_sp()
decov_obj.deconvolution(
    sketch_dim=512,           # randomized-sketch dimension; raise to 1024 for HD-2µm
    lambda_spatial=5000,      # spatial smoothness; raise to 10000 for sparse / noisy data
    n_hvg=2000,               # raise to 3000 when accuracy matters more than speed
)

# Results — three accessors (same payload, different wrappings):
decov_obj.adata_cell2location               # AnnData with proportions in .X
decov_obj.adata_sp.obsm['flashdeconv']      # DataFrame of proportions
decov_obj.adata_sp.obs['flashdeconv_dominant']  # dominant cell type per spot

# Spatial visualisation reuses the same plot stack as cell2location
ov.pl.plot_spatial(
    adata=decov_obj.adata_cell2location,
    color=clust_labels, labels=clust_labels,
    show_img=True, style='fast',
    max_color_quantile=0.992,
    circle_diameter=4,
    colorbar_position='right',
    palette=color_dict,
)
```

When to pick FlashDeconv over the others:
| Feature | FlashDeconv | Tangram | cell2location |
|---|---|---|---|
| GPU required | No | Optional | Recommended |
| 10k spots wall-time | ~2 min | ~15 min | ~60 min |
| Native Visium HD support | Yes | Limited | Limited |
| Built-in spatial regularisation | Yes (`lambda_spatial`) | No | No |

Parameter tuning quick reference (from the FlashDeconv tutorial):
- Noisy / sparse data → raise `lambda_spatial` to 10000.
- Visium HD 2 µm → raise `sketch_dim` to 1024.
- Accuracy-critical → raise `n_hvg` to 3000.

### Starfysh Archetypal Deconvolution

Treat the file names below as placeholders for your own local exports; the matching public spatial tutorials describe the expected AnnData counts matrix and signature table formats.

```python
from omicverse.external.starfysh import AA, utils, plot_utils

visium_args = utils.prepare_data(adata_path="data/counts.h5ad",
                                  signature_path="data/signatures.csv",
                                  min_cells=10, filter_hvg=True, n_top_genes=3000)
adata, adata_normed = visium_args.get_adata()
aa_model = AA.ArchetypalAnalysis(adata_orig=adata_normed)
aa_model.fit(k=12, n_init=10)
visium_args = utils.refine_anchors(visium_args, aa_model, add_marker=True)
model, history = utils.run_starfysh(visium_args, poe=False, n_repeat=5, lr=5e-3, max_epochs=500)
```

## Stage 3: Downstream Analysis

### Spatial Clustering

```python
# GraphST + mclust
ov.utils.cluster(adata, use_rep='graphst|original|X_pca', method='mclust', n_components=7)

# STAGATE
ov.utils.cluster(adata, use_rep='STAGATE', method='mclust', n_components=7)

# Merge small clusters
ov.space.merge_cluster(adata, groupby='mclust', resolution=0.5)
```

Algorithm choice: GraphST and STAGATE require precalculated latent spaces. For standard clustering without spatial-aware embeddings, use Leiden/Louvain directly.

### Multi-Slice Integration (STAligner)

```python
import anndata as ad
Batch_list = [ov.read(p) for p in slice_paths]
adata_concat = ad.concat(Batch_list, label='slice_name', keys=section_ids)
STAligner_obj = ov.space.pySTAligner(adata=adata_concat, batch_key='slice_name',
                                      hidden_dims=[256, 64], use_gpu=True)
STAligner_obj.train_STAligner_subgraph(nepochs=800, lr=1e-3)
STAligner_obj.train()
adata_aligned = STAligner_obj.predicted()
sc.pp.neighbors(adata_aligned, use_rep='STAligner')
```

### Spatial Trajectories

**SpaceFlow** — pseudo-spatial maps:
```python
sf_obj = ov.space.pySpaceFlow(adata)
sf_obj.train(spatial_regularization_strength=0.1, num_epochs=300, patience=50)
sf_obj.cal_pSM(n_neighbors=20, resolution=1.0)
```

**STT** — transition dynamics:
```python
STT_obj = ov.space.STT(adata, spatial_loc='xy_loc', region='Region', n_neighbors=20)
STT_obj.stage_estimate()
STT_obj.train(n_states=9, n_iter=15, weight_connectivities=0.5)
```

### Cell Communication (COMMOT + FlowSig)

```python
df_cellchat = ov.external.commot.pp.ligand_receptor_database(species='human', database='cellchat')
df_cellchat = ov.external.commot.pp.filter_lr_database(df_cellchat, adata, min_expr_frac=0.05)
ov.external.commot.tl.spatial_communication(adata, lr_database=df_cellchat,
                                              distance_threshold=500, result_prefix='cellchat')
# FlowSig network
adata.layers['normalized'] = adata.X.copy()
ov.external.flowsig.tl.construct_intercellular_flow_network(
    adata, commot_output_key='commot-cellchat',
    flowsig_output_key='flowsig-cellchat', edge_threshold=0.7)
```

### Structural Layers (GASTON) and Slice Alignment (SLAT)

**GASTON** — iso-depth estimation:
```python
gas_obj = ov.space.GASTON(adata)
A = gas_obj.prepare_inputs(n_pcs=50)
gas_obj.load_rescale(A)
gas_obj.train(hidden_dims=[64, 32], dropout=0.1, max_epochs=2000)
gaston_isodepth, gaston_labels = gas_obj.cal_iso_depth(n_layers=5)
```

**SLAT** — cross-slice alignment:
```python
from omicverse.external.scSLAT.model import Cal_Spatial_Net, load_anndatas, run_SLAT, spatial_match
Cal_Spatial_Net(adata1, k_cutoff=20, model='KNN')
Cal_Spatial_Net(adata2, k_cutoff=20, model='KNN')
edges, features = load_anndatas([adata1, adata2], feature='DPCA', check_order=False)
embeddings, *_ = run_SLAT(features, edges, LGCN_layer=5)
best, index, distance = spatial_match(embeddings, adatas=[adata1, adata2])
```

## Troubleshooting

- **`ValueError: spatial coordinates out of bounds` after rotation**: Cast `adata.obsm['spatial']` to `float64` BEFORE calling `rotate_space_visium`. Integer coordinates lose precision during trigonometric rotation.
- **Cellpose segmentation fails with memory error**: For large `.btf` images, use `backend='tifffile'` to memory-map the image. Reduce `buffer` parameter if GPU memory is insufficient.
- **Gene ID overlap failure in Tangram/cell2location**: Harmonise identifiers (ENSEMBL vs gene symbols) between `adata_sc` and `adata_sp` before calling `preprocess_sc`/`preprocess_sp`. Drop non-overlapping genes.
- **`mclust` clustering error**: Requires `rpy2` and the R `mclust` package. If R bindings are unavailable, switch to `method='louvain'` or `method='leiden'`.
- **STAligner/SpaceFlow embeddings collapse to a single point**: Verify `adata.obsm['spatial']` exists and coordinates are scaled appropriately. Tune learning rate (try `lr=5e-4`) and regularisation strength.
- **FlowSig returns empty network**: Build spatial neighbor graphs before Moran's I filtering. Increase bootstraps or lower `edge_threshold` (try 0.5) if the network is too sparse.
- **GASTON `RuntimeError` in training**: Provide a writable `out_dir` path. PyTorch nondeterminism may cause variation between runs—set `torch.manual_seed()` for reproducibility.
- **SLAT alignment has many low-quality matches**: Regenerate spatial graphs with a higher `k_cutoff` value. Inspect `low_quality_index` flags and filter cells with high distance scores.
- **STT pathway enrichment fails**: `gseapy` needs network access for gene set downloads. Cache gene sets locally with `ov.utils.geneset_prepare()` and pass the dictionary directly.

## Dependencies
- Core: `omicverse`, `scanpy`, `anndata`, `squidpy`, `numpy`, `matplotlib`
- Segmentation: `cellpose`, `opencv-python`/`tifffile`, optional GPU PyTorch
- Deconvolution: `tangram-sc`, `cell2location`, `pytorch-lightning`; Starfysh needs `torch`, `scikit-learn`
- Downstream: `scikit-learn`, `commot`, `flowsig`, `gseapy`, torch-backed modules (STAligner, SpaceFlow, GASTON, SLAT)

## Examples
- "Crop and rotate my Visium slide, then run cellpose segmentation on the HD data and aggregate to cell-level AnnData."
- "Deconvolve my lymph node spatial data with Tangram and cell2location, compare proportions, and plot cell-type maps."
- "Integrate three DLPFC slices with STAligner, cluster with STAGATE, and infer communication with COMMOT+FlowSig."

## References
- Quick copy/paste commands: [`reference.md`](reference.md)
