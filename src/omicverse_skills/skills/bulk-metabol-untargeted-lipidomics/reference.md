# Untargeted LC-MS & Lipidomics — quick commands

## LC-MS ingest + mummichog

```python
import omicverse as ov
import numpy as np
ov.plot_set()

# 1) Load LC-MS peak table — m/z + RT parsed into adata.var
adata = ov.metabol.read_lcms(
    'malaria_feature_table.csv',
    feature_id_sep='__',          # whatever your file uses; '/' is common
    label_row='Label',
    transpose=True,
)
# adata.var['m_z'] (Da), adata.var['rt'] (min); obs['group'] from label_row

# 2) Preprocess — DON'T impute LC-MS zeros (they're below-detection)
adata = ov.metabol.normalize(adata, method='pqn')
adata = ov.metabol.transform(adata, method='log')

# 3) Differential — per-peak p-values feed mummichog
deg = ov.metabol.differential(
    adata, group_col='group',
    group_a='Semi_immue', group_b='Naive',
    method='welch_t', log_transformed=True,
)

# 4) Honest volcano on raw p-values (BH-FDR is too strict on untargeted LC-MS)
fig, ax = ov.metabol.volcano(
    deg, padj_thresh=0.01, log2fc_thresh=2.0,
    label_top_n=0, use_pvalue=True, clip_log2fc=8.0,
)

# 5) Pre-fetch DBs once to avoid re-network if you call mummichog multiple times
pathways = ov.metabol.load_pathways()           # ~550 KEGG pathways, cached
mass_db  = ov.metabol.fetch_chebi_compounds()   # ~54k ChEBI compounds, cached

# 6) Mummichog — pure-Python port; permutation-based
mumm = ov.metabol.mummichog_basic(
    mz=adata.var['m_z'].values,
    pvalue=deg['pvalue'].values,        # raw, not padj
    polarity='positive',
    ppm=10.0,
    significance_cutoff=0.05,
    n_perm=1000, min_overlap=2,
    mass_db=mass_db, pathways=pathways,
    seed=0,
)
print(mumm.head(10)[['pathway', 'overlap', 'set_size', 'pvalue', 'padj']])

# 7) Optional — peak-level adduct annotation
ann = ov.metabol.annotate_peaks(
    adata.var['m_z'].values,
    polarity='positive', ppm=10.0,
    mass_db=mass_db,
)
print(f'{ann["peak_idx"].nunique()} peaks have ≥1 KEGG candidate at 10 ppm')
```

## Mummichog sanity check (seed a known pathway)

```python
# Build a synthetic m/z set from a known KEGG pathway and confirm recovery.
pathways = ov.metabol.load_pathways()
mass_db  = ov.metabol.fetch_chebi_compounds()

tca_ids  = pathways['Citrate cycle (TCA cycle)']
tca = mass_db[mass_db['mw'].notna() & mass_db['kegg'].isin(tca_ids)].reset_index(drop=True)

rng = np.random.default_rng(0)
hit_mz = tca['mw'].to_numpy() + 1.00728            # synthesize [M+H]+
bg_mz  = rng.uniform(50, 1200, size=80)
all_mz = np.concatenate([hit_mz, bg_mz])
pvals  = np.concatenate([np.full(len(hit_mz), 0.001),
                         np.full(80, 0.5)])

sanity = ov.metabol.mummichog_basic(
    mz=all_mz, pvalue=pvals, polarity='positive',
    ppm=10.0, significance_cutoff=0.05, n_perm=500, min_overlap=2,
    mass_db=mass_db, pathways=pathways,
)
print(sanity.head(5)[['pathway', 'overlap', 'set_size', 'pvalue', 'padj']])
# 'Citrate cycle (TCA cycle)' should appear with low p-value
```

## External mummichog (full Li reference implementation)

```python
# pip install mummichog
mumm = ov.metabol.mummichog_external(
    mz=adata.var['m_z'].values,
    pvalue=deg['pvalue'].values,
    polarity='positive', ppm=10.0,
    significance_cutoff=0.05,
)
# Adds activity-network scoring + LIBSDB / SMPDB pathways out of the box.
```

## Lipidomics — parse + annotate + class aggregation

```python
import omicverse as ov
import pandas as pd, anndata as ad

# Build AnnData
matrix = pd.read_csv('brca_matrix.csv', index_col=0)
clin   = pd.read_csv('brca_clin.csv',  index_col=0)
adata = ad.AnnData(matrix.T)
adata.obs = clin.loc[adata.obs_names].copy()
adata.obs['group'] = adata.obs['Group']

# Parse a few names to confirm the regex understands your format
for name in ['PC 34:1', 'TAG 54:3', 'Cer d18:1/24:0', 'LPC 18:0', 'TAG 54:3;O']:
    print(name, '→', ov.metabol.parse_lipid(name))

# Annotate every var; unparseable rows get NaN in lipid_class
adata = ov.metabol.annotate_lipids(adata)
print(adata.var['lipid_class'].value_counts())

# Class-level matrix (optional)
class_adata = ov.metabol.aggregate_by_class(adata, agg='sum')
class_adata.obs = adata.obs        # propagate clinical metadata
```

## Lipidomics — species-level differential + LION enrichment

```python
# Drop unparseable species before stats
adata_clean = adata[:, adata.var['lipid_class'].notna()].copy()

deg = ov.metabol.differential(
    adata_clean, group_col='group',
    group_a='Cancer', group_b='Benign',
    method='welch_t', log_transformed=True,
)
fig, ax = ov.metabol.volcano(deg, padj_thresh=0.10, log2fc_thresh=0.3, label_top_n=8)

# LION ontology enrichment — feeds species-level NAMES, not classes
hits = deg[deg.padj < 0.10].index.tolist()
background = deg.index.tolist()
lion = ov.metabol.lion_enrichment(hits, background, min_size=2)
print(lion.head(10)[['term', 'category', 'overlap', 'set_size',
                     'odds_ratio', 'pvalue', 'padj']])
```

## Pre-fetch the LION ontology when running multiple enrichments

```python
ontology = ov.metabol.fetch_lion_associations()
lion = ov.metabol.lion_enrichment(hits, background, ontology=ontology, min_size=3)
```

## LipidIdentity convenience checks

```python
lid = ov.metabol.parse_lipid('PUFA-rich PC 38:6')   # not real LIPID MAPS, just shape
lid = ov.metabol.parse_lipid('PC 38:6')
print(lid.lipid_class, lid.total_carbons, lid.total_db)
print('saturated?         ', lid.is_saturated())
print('polyunsaturated?   ', lid.is_polyunsaturated(threshold=2))
```
