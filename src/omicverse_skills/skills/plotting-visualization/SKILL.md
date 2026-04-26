---
name: omicverse-visualization-for-bulk-color-systems-and-single-cell-d
title: OmicVerse visualization for bulk, color systems, and single-cell data
description: "OmicVerse plotting: volcano, venn, boxplot, embedding, density, dotplot, convex hull, stacked bar, and Forbidden City color palettes."
---

# OmicVerse visualization for bulk, color systems, and single-cell data

## Overview
Leverage this skill when a user wants help recreating or adapting plots from the OmicVerse plotting tutorials:
- [`t_visualize_bulk.ipynb`](https://omicverse.readthedocs.io/en/latest/Tutorials-plotting/t_visualize_bulk/)
- [`t_visualize_colorsystem.ipynb`](https://omicverse.readthedocs.io/en/latest/Tutorials-plotting/t_visualize_colorsystem/)
- [`t_visualize_single.ipynb`](https://omicverse.readthedocs.io/en/latest/Tutorials-plotting/t_visualize_single/)

It covers how to configure OmicVerse's plotting style, choose colors from the Forbidden City palette, and generate bulk as well as single-cell specific figures.

## Instructions
1. **Set up the plotting environment**
   - Import `omicverse as ov`, `matplotlib.pyplot as plt`, and other libraries required by the user's request (`pandas`, `seaborn`, `scanpy`, etc.).
   - Call `ov.ov_plot_set()` (or `ov.plot_set()` depending on the installed version) to apply OmicVerse's default styling before generating figures.
   - Load example data via `ov.read(...)`/`ov.pp.preprocess(...)` or instruct users to supply their own AnnData/CSV files.
2. **Bulk RNA-seq visuals (`t_visualize_bulk`)**
   - Use `ov.pl.venn(sets=..., palette=...)` to display overlaps among DEG lists (no more than 4 groups). Encourage setting `sets` as a dictionary of set names → gene lists.
   - For volcano plots, load the DEG table (`result = ov.read('...csv')`) and call `ov.pl.volcano(result, pval_name='qvalue', fc_name='log2FoldChange', ...)`. Explain optional keyword arguments such as `sig_pvalue`, `sig_fc`, `palette`, and label formatting.
   - To compare group distributions with box plots, gather long-form data (e.g., from `seaborn.load_dataset('tips')`) and invoke `ov.pl.boxplot(data, x_value=..., y_value=..., hue=..., ax=ax, palette=...)`. Mention how to adjust figure size, legend placement, and significance annotations.
3. **Color management (`t_visualize_colorsystem`)**
   - Introduce the color book via `fb = ov.pl.ForbiddenCity()` and demonstrate `fb.get_color(name='凝夜紫')` for specific hues.
   - Show how to pull predefined palettes (`ov.pl.green_color`, `ov.pl.red_color`, etc.) and build dicts mapping cell types/groups to color hex codes.
   - For segmented gradients, combine colors and call `ov.pl.get_cmap_seg(colors, name='custom')`, then pass the colormap into Matplotlib/Scanpy plotting functions.
   - Highlight using these palettes in embeddings: `ov.pl.embedding(adata, basis='X_umap', color='clusters', palette=color_dict, ax=ax)`.
4. **Single-cell visualizations (`t_visualize_single`)**
   - Remind users to preprocess AnnData if needed (`adata = ov.pp.preprocess(adata, mode='shiftlog|pearson', n_HVGs=2000)`).
   - **IMPORTANT - Data validation**: Before plotting, always verify that required data exists:
     ```python
     # Before plotting by clustering or other categorical variable
     color_col = 'leiden'  # or 'clusters', 'celltype', etc.
     if color_col not in adata.obs.columns:
         raise ValueError(f"Column '{color_col}' not found in adata.obs. Available columns: {list(adata.obs.columns)}")

     # Before plotting embeddings
     basis = 'X_umap'  # or 'X_pca', 'X_tsne', etc.
     if basis not in adata.obsm.keys():
         raise ValueError(f"Embedding '{basis}' not found in adata.obsm. Available embeddings: {list(adata.obsm.keys())}")
     ```
   - For palette optimization, use `ov.pl.optim_palette(adata, basis='X_umap', colors='clusters')` to auto-generate color schemes when categories clash.
   - Reproduce stacked proportions with `ov.pl.cellproportion(adata, groupby='clusters', celltype_clusters='celltype', ax=ax)` and transform into stacked area charts by setting `kind='area'`.
   - Showcase compound embedding utilities:
     - `ov.pl.embedding_celltype` to place counts/proportions alongside UMAPs.
     - `ov.pl.ConvexHull` or `ov.pl.contour` for highlighting regions of interest.
     - `ov.pl.embedding_adjust` to reposition legends automatically.
     - `ov.pl.embedding_density` for density overlays, controlling smoothness with `adjust`.
   - For spatial gene density, describe the workflow: `ov.pl.calculate_gene_density(adata, genes=[...], basis='spatial')`, then overlay with `ov.pl.embedding(..., layer='gene_density', cmap='...')`.
   - Cover additional charts like `ov.pl.single_group_boxplot`, `ov.pl.bardotplot`, `ov.pl.dotplot`, and `ov.pl.marker_heatmap`, emphasizing input formats (long-form DataFrame vs. AnnData with `.obs` annotations) and optional helpers such as `ov.pl.add_palue` for manual p-value annotations.
   - **Circular plot1cell view (`ov.pl.plot1cell`)**: a Wu 2021 plot1cell-style circular UMAP with concentric metadata tracks. Clusters become arc sectors on the circumference (sector length ∝ `log10(n_cells)`); the embedding is drawn inside the unit circle with a Gaussian-KDE contour overlay; each entry in `tracks` becomes one extra ring coloured by the run-length segments of that metadata column within each cluster sector. Single call:
     ```python
     ov.pl.plot1cell(
         adata, clusters='cell_type', basis='X_umap',
         tracks=['compartment', 'tissue', 'sex', 'age'],
         point_size=2, point_alpha=0.35,
         figsize=(10, 10), label_fontsize=7,
     )
     ```
     - **Scale `point_size` / `point_alpha` to cohort size**: ~10k cells `point_size=6, alpha=0.5`; 50k `2, 0.35`; 100k `1, 0.25`; 200k+ `0.8, 0.2`.
     - **`tracks` order matters**: each track becomes a concentric ring; the first listed sits closest to the centre. Put the most-relevant metadata first.
     - **`basis`**: pass any 2-D embedding key (`X_umap`, `X_tSNE`, `X_umap_Harmony_scDonor_snBatch`, etc.) — the function accepts arbitrary obsm keys.
     - **Memory tip**: plot1cell only needs `obsm[basis]` + `obs[clusters]` + `obs[tracks]`; for 200k+ cohorts, slice `adata = adata[:, :200].copy()` to a few hundred genes before plotting (the gene matrix is irrelevant to this view).
     - **Custom palettes**: pass `cluster_palette` (cluster-arc colours) and/or `track_palettes` (list of per-track palettes) to override defaults.
5. **Finishing touches and exports**
   - Encourage adding titles, axis labels, and `fig.tight_layout()` to prevent clipping.
   - Suggest saving figures with `fig.savefig('plot.png', dpi=300, bbox_inches='tight')` and documenting color mappings for reproducibility.
   - Troubleshoot common issues:
     - **Missing AnnData keys**: Always validate `adata.obs` columns and `adata.obsm` embeddings exist before plotting
     - **Palette names not found**: Verify color dictionaries match actual category values
     - **Matplotlib font rendering**: When using Chinese characters, ensure appropriate fonts are installed
     - **"Could not find X in adata.obs"**: Check that clustering or annotation has been performed before trying to visualize results. Use defensive checks to compute missing prerequisites on-the-fly.

## Examples
- "Plot a three-set Venn diagram of overlapping DEG lists and reuse Forbidden City colors for consistency."
- "Load the dentate gyrus AnnData, color clusters with `fb.get_color` selections, and render an embedding with adjusted legend placement."
- "Generate single-cell proportion bar/area plots plus gene-density overlays using OmicVerse helper functions."
- "Build a `plot1cell` circular view of a 100k-cell cohort with disease / tissue / assay / sex tracks and `point_size=1, point_alpha=0.25`."

## References
- Bulk tutorial: [`t_visualize_bulk.ipynb`](https://omicverse.readthedocs.io/en/latest/Tutorials-plotting/t_visualize_bulk/)
- Color system tutorial: [`t_visualize_colorsystem.ipynb`](https://omicverse.readthedocs.io/en/latest/Tutorials-plotting/t_visualize_colorsystem/)
- Single-cell tutorial: [`t_visualize_single.ipynb`](https://omicverse.readthedocs.io/en/latest/Tutorials-plotting/t_visualize_single/)
- plot1cell circular view: [`t_plot1cell.ipynb`](https://omicverse.readthedocs.io/en/latest/Tutorials-plotting/t_plot1cell/)
- Quick reference snippets: [`reference.md`](reference.md)
