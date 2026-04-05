# Source Grounding

## Inspected Interfaces

- `ov.single.SCENIC(adata, db_glob, motif_path, n_jobs=8)`
- `SCENIC.cal_grn(method='regdiffusion'|'grnboost2'|'genie3', layer='counts', tf_names=None, **kwargs)`
- `SCENIC.cal_regulons(rho_mask_dropouts=True, seed=42, **kwargs)`
- `modules_from_adjacencies(adjacencies, ex_mtx, thresholds=(0.75, 0.90), top_n_targets=(50,), top_n_regulators=(5, 10, 50), min_genes=20, absolute_thresholds=False, rho_dichotomize=True, keep_only_activating=True, rho_threshold=0.03, rho_mask_dropouts=False)`
- `regulon_specificity_scores(auc_mtx, cell_type_series)`
- `binarize(auc_mtx, threshold_overides=None, seed=None, num_workers=1)`
- `plot_rss(rss, cell_type, top_n=5, max_n=None, ax=None)`
- `sc.tl.rank_genes_groups(..., method='logreg'|'t-test'|'wilcoxon'|'t-test_overestim_var', use_raw=..., layer=...)`
- `ov.single.cosg(..., mu=1, remove_lowly_expressed=False, expressed_pct=0.1, key_added=None, calculate_logfoldchanges=True, use_raw=True, layer=None)`
- `ov.single.build_correlation_network_umap_layout(embedding_df, correlation_threshold=0.7, umap_neighbors=15, min_dist=0.1)`
- `ov.single.add_tf_regulation(G, tf_gene_dict)`
- `ov.single.plot_grn(G, pos, tf_list, temporal_df, tf_gene_dict, figsize=(6, 6), top_tf_target_num=5, title='GRN', ax=None, fontsize=12, cmap='RdBu_r')`

## Source Notes

- The live `SCENIC` constructor checks that the ranking databases and motif annotation table exist before any inference starts.
- `SCENIC.cal_grn(...)` raises if the requested layer is missing or already appears log-normalized.
- The `method` branch in `SCENIC.cal_grn(...)` currently supports `regdiffusion`, `grnboost2`, and `genie3`.
- `SCENIC.cal_regulons(...)` does not expose every pruning parameter directly in its signature; it forwards extra kwargs to pySCENIC module construction.
- pySCENIC `modules_from_adjacencies(...)` combines percentile-threshold modules, top-target modules, and top-regulator modules before pruning.
- The live pySCENIC default for `rho_mask_dropouts` is `False`, but the notebook sets it to `True`.
- `binarize(...)` uses the parameter name `threshold_overides` in the installed source.
- `plot_rss(...)` annotates only one cell type per call and labels the top regulons selected by `top_n`.
- `scanpy.tl.rank_genes_groups(...)` supports `logreg`, `t-test`, `wilcoxon`, and `t-test_overestim_var`, and it disallows `layer` together with `use_raw=True`.
- `ov.single.cosg(...)` exposes its own branch controls through `mu`, `remove_lowly_expressed`, `expressed_pct`, `calculate_logfoldchanges`, `use_raw`, and `layer`.

## Notebook-Grounded Claims

- The notebook uses a tutorial-era raw-count layer name instead of the live `counts` default.
- The notebook's regulon branch passes `rho_mask_dropouts=True`, `thresholds=(0.75, 0.9)`, `top_n_targets=(50,)`, and `top_n_regulators=(5, 10, 50)`.
- The notebook's downstream interpretation stages are RSS, regulon embedding overlays, regulon-marker ranking, binarization, and TF-centered GRN plotting.
