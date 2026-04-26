# Pathway / multifactor / DGCA / MOFA / MTBLS1 — quick commands

## MSEA pathway enrichment (ORA + GSEA)

```python
import omicverse as ov
import matplotlib.pyplot as plt
ov.plot_set()

# preprocessed (PQN+log) AnnData
adata = ov.metabol.read_metaboanalyst('human_cachexia.csv', group_col='Muscle loss')
adata = ov.metabol.normalize(adata, method='pqn')
adata = ov.metabol.transform(adata, method='log')

deg = ov.metabol.differential(adata, group_col='group',
                              group_a='cachexic', group_b='control',
                              method='welch_t', log_transformed=True)

# Coverage check (optional)
ids = ov.metabol.map_ids(adata.var_names.tolist(),
                         targets=('hmdb', 'kegg', 'chebi'))

# ORA — needs a hit list + background
hits = deg[deg.padj < 0.20].index.tolist()
background = deg.index.tolist()
ora = ov.metabol.msea_ora(hits, background, min_size=3)
ov.metabol.pathway_bar(ora, score_col='pvalue', top_n=10)
ov.metabol.pathway_dot(ora, size_col='overlap', x_col='odds_ratio',
                       color_col='pvalue', top_n=10)

# GSEA — pathway names in the INDEX
gsea = ov.metabol.msea_gsea(deg, stat_col='stat', n_perm=500, seed=0)
gsea_plot = gsea.reset_index().rename(columns={
    gsea.index.name or 'index': 'pathway',
})
ov.metabol.pathway_dot(gsea_plot, size_col='matched_size', x_col='nes',
                       color_col='pval', top_n=10)
plt.show()
```

## VIP × DEG shortlist

```python
adata_pareto = ov.metabol.transform(adata, method='pareto', stash_raw=False)
opls = ov.metabol.opls_da(adata_pareto, group_col='group', n_ortho=1)
vip = opls.to_vip_table(adata.var_names)

combined = deg.join(vip[['vip']]).sort_values('vip', ascending=False)
shortlist = combined[(combined['padj'] < 0.20) & (combined['vip'] > 1.0)]
print(shortlist[['pvalue', 'padj', 'log2fc', 'vip']])
```

## Multi-factor — ASCA

```python
asca_res = ov.metabol.asca(
    adata, factors=['treatment', 'time'],
    include_interactions=True,
    n_components=2, n_permutations=500, seed=0,
)
print(asca_res.summary())
ov.metabol.asca_variance_bar(asca_res); plt.show()
```

## Multi-factor — mixed model (repeated measures)

```python
tbl = ov.metabol.mixed_model(
    adata,
    formula='treatment + time',
    groups='patient',
    term='treatment[T.drug]',     # categorical contrast — use exact patsy syntax
)
print(tbl.sort_values('pvalue').head(10))
```

## Multi-factor — MEBA two-group time-series

```python
# obs columns required: 'group', 'time', 'subject'
meba_tbl = ov.metabol.meba(adata_ts, group_col='group',
                            time_col='time', subject_col='subject')
print(meba_tbl.sort_values('pvalue').head(10))
```

## DGCA + per-condition correlation networks

```python
# differential correlation between groups
dc = ov.metabol.dgca(
    adata, group_col='group',
    group_a='cachexic', group_b='control',
    method='spearman',
    abs_r_threshold=0.3,
)
ov.metabol.dgca_class_bar(dc); plt.show()

# per-condition networks — keep both for set arithmetic
edges_a = ov.metabol.corr_network(adata, group_col='group', group='cachexic',
                                  method='spearman',
                                  abs_r_threshold=0.5, padj_threshold=0.05)
edges_b = ov.metabol.corr_network(adata, group_col='group', group='control',
                                  method='spearman',
                                  abs_r_threshold=0.5, padj_threshold=0.05)
print(f'cachexic n_samples: {edges_a.attrs["n_samples"]}, edges: {len(edges_a)}')
print(f'control  n_samples: {edges_b.attrs["n_samples"]}, edges: {len(edges_b)}')

# shared / specific edges
def key(row): return frozenset((row['feature_a'], row['feature_b']))
keys_a = set(edges_a.apply(key, axis=1))
keys_b = set(edges_b.apply(key, axis=1))
print(f'shared: {len(keys_a & keys_b)}  A-only: {len(keys_a - keys_b)}  B-only: {len(keys_b - keys_a)}')

# plot top edges
ov.metabol.corr_network_plot(edges_a.head(40),
                             figsize=(7, 6), label_font_size=6)
plt.show()
```

## Multi-omics MOFA+

```python
# Both views must share obs_names in the same order
factors = ov.metabol.run_mofa(
    views={'metabol': adata_m, 'rna': adata_r},
    n_factors=5,
    outfile='mofa_demo.hdf5',
    scale_views=True,
    max_iter=200,
    convergence_mode='fast',     # 'medium'/'slow' for reportable
    seed=0,
)
print(f'retained {factors.shape[1]} factors')
print(factors.head())
```

## MTBLS1 worked example (T2D urine NMR)

```python
import omicverse as ov
import matplotlib.pyplot as plt
ov.plot_set()

# 1) Ingest from MetaboLights public mirror
adata = ov.utils.load_metabolights(
    'MTBLS1',
    group_col='Factor Value[Metabolic syndrome]',
    cache_dir='mtbls1_demo',
)
print('group counts:', adata.obs['group'].value_counts().to_dict())

# 2) Sample QC
qc = ov.metabol.sample_qc(adata, n_components=3, alpha=0.95)
print(f'flagged: {qc["is_outlier"].sum()} / {len(qc)}')
ov.metabol.sample_qc_plot(qc); plt.show()

# 3) Preprocess: CV-filter → impute → normalize → transform
adata = ov.metabol.cv_filter(adata, across='all', cv_threshold=1.5)
adata = ov.metabol.impute(adata, method='qrilc', seed=0)
adata = ov.metabol.normalize(adata, method='pqn')
adata_log = ov.metabol.transform(adata, method='log')
adata_pareto = ov.metabol.transform(adata_log, method='pareto', stash_raw=False)

# 4) Univariate + volcano (log)
deg = ov.metabol.differential(adata_log, group_col='group',
                              group_a='diabetes mellitus', group_b='Control Group',
                              method='welch_t', log_transformed=True)
ov.metabol.volcano(deg, padj_thresh=0.05, log2fc_thresh=0.5, label_top_n=8); plt.show()

# 5) Multivariate (Pareto) + S-plot + VIP with named labels
opls = ov.metabol.opls_da(adata_pareto, group_col='group',
                          group_a='diabetes mellitus', group_b='Control Group',
                          n_ortho=1)
ov.metabol.s_plot(opls, adata_pareto, label_top_n=8); plt.show()
mb_names = adata_pareto.var['metabolite_identification'].values
ov.metabol.vip_bar(opls, mb_names, top_n=15); plt.show()

# 6) MSEA ORA on names
mass_db = ov.metabol.fetch_chebi_compounds()
hit_idx = deg[deg['padj'] < 0.05].index
hit_names = adata.var.loc[hit_idx, 'metabolite_identification'].tolist()
bg_names = adata.var['metabolite_identification'].dropna().tolist()
enr = ov.metabol.msea_ora(hit_names, bg_names, mass_db=mass_db)
ov.metabol.pathway_bar(enr, top_n=10); plt.show()

# 7) ASCA over phenotype + sex
asca_res = ov.metabol.asca(
    adata_pareto, factors=['group', 'Factor Value[Gender]'],
    include_interactions=True, n_permutations=100, seed=0,
)
print(asca_res.summary().round(4))
ov.metabol.asca_variance_bar(asca_res); plt.show()

# 8) ROC + biomarker panel (delegated to multivariate skill)
auc = ov.metabol.roc_feature(adata_log, group_col='group',
                             pos_group='diabetes mellitus', neg_group='Control Group')
panel = ov.metabol.biomarker_panel(
    adata_log, group_col='group',
    pos_group='diabetes mellitus', neg_group='Control Group',
    features=10, classifier='lr',
    cv_outer=5, cv_inner=3, n_permutations=100, seed=0,
)

# 9) DGCA + named correlation network
dc = ov.metabol.dgca(adata_log, group_col='group',
                     group_a='diabetes mellitus', group_b='Control Group',
                     method='spearman')
ov.metabol.dgca_class_bar(dc); plt.show()

edges_t2d = ov.metabol.corr_network(
    adata_log, group_col='group', group='diabetes mellitus',
    method='spearman', abs_r_threshold=0.5, padj_threshold=0.05,
)
named = adata.var['metabolite_identification'].to_dict()
edges_plot = edges_t2d.copy()
edges_plot['feature_a'] = edges_plot['feature_a'].map(named).fillna(edges_plot['feature_a'])
edges_plot['feature_b'] = edges_plot['feature_b'].map(named).fillna(edges_plot['feature_b'])
ov.metabol.corr_network_plot(edges_plot.head(40), figsize=(7, 6), label_font_size=6)
plt.show()
```
