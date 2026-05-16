---
name: omicverse-bulk-metabol-untargeted-lipidomics
description: Two adjacent LC-MS workflows on AnnData — (1) untargeted metabolomics with m/z-based peak annotation, mummichog pathway inference and adduct-ppm matching, and (2) lipidomics with LIPID MAPS shorthand parsing, lipid-class aggregation, and LION term enrichment. Use when converting `t_metabol_04_untargeted` or `t_metabol_05_lipidomics` into a reusable skill, when the input feature IDs encode `m/z`/`RT`, or when the var_names look like `PC 34:1` / `Cer d18:1/24:0` / `TAG 54:3`.
---

# OmicVerse Bulk Metabolomics — Untargeted LC-MS & Lipidomics

## Goal

Take an LC-MS peak table (m/z + RT-coded features) or a lipidomics matrix (LIPID MAPS shorthand var-names) into an `AnnData`, then run the two ID-aware downstream paths — peak annotation + mummichog pathway inference for untargeted metabolomics, lipid-class parsing + class aggregation + LION enrichment for lipidomics. The output is a list of pathway / ontology hits backed by the empirical feature p-values, plus a class-aggregated lipid matrix when requested.

This skill assumes the upstream preprocessing chain (impute / normalize / log) has already run; for the chain see the `omicverse-bulk-metabol-preprocessing` skill. ID mapping for *named* metabolites (HMDB / KEGG / ChEBI ID resolution) and MSEA pathway enrichment are in the pathway-multifactor skill.

## Quick Workflow

### Untargeted LC-MS path (t_metabol_04)

1. Load the LC-MS peak table with `ov.metabol.read_lcms(..., feature_id_sep=...)` so `m/z` and `RT` get parsed into `var['m_z']` and `var['rt']`.
2. Apply preprocessing (PQN + log) — **don't impute zeros** for LC-MS (zeros = below-detection).
3. Run `differential` to get a per-peak p-value vector.
4. Untargeted LC-MS often has thousands of peaks vs. low n — BH-FDR is too strict to discriminate. For the volcano, pass `use_pvalue=True` and clip the `log2fc` axis to keep below-detection extremes from blowing the figure.
5. Enrich pathways from m/z + p-value with `ov.metabol.mummichog_basic(...)` — a permutation-based pathway test that handles the FDR burden internally.
6. Optionally inspect candidate compound annotations per peak with `ov.metabol.annotate_peaks(...)` — adduct-resolved KEGG candidates by ppm tolerance.
7. For sanity checks: build a synthetic m/z set from a known KEGG pathway (e.g. via `fetch_chebi_compounds()` + `load_pathways()`), seed it with low p-values, and confirm mummichog recovers the seeded pathway.

### Lipidomics path (t_metabol_05)

1. Build the `AnnData` directly from a wide CSV — `ov.metabol.read_wide(...)` or `anndata.AnnData(df)` then attach clinical metadata.
2. Sanity-check the parser on a few var_names: `ov.metabol.parse_lipid('PC 34:1')` returns a `LipidIdentity` (`lipid_class`, `total_carbons`, `total_db`, `backbone`).
3. Run `ov.metabol.annotate_lipids(adata)` — adds `var['lipid_class']`, `var['total_carbons']`, `var['total_db']`. Unparseable rows get `NaN`; drop them before downstream stats.
4. (Optional) `ov.metabol.aggregate_by_class(adata, agg='sum')` for a class-level matrix (`n_vars = n_lipid_classes`) — useful for class-resolved QC and parsimonious modelling.
5. Run `differential` on the cleaned (non-NaN) species-level matrix to get per-lipid hits.
6. Enrich LION lipid-ontology terms from the hit list with `ov.metabol.lion_enrichment(hits, background, min_size=2)`. Term categories include lipid class, function, fatty-acid composition.

## Interface Summary

LC-MS ingest (creates `AnnData` with `var['m_z']`, `var['rt']`):
- `ov.metabol.read_lcms(path, *, feature_id_sep='/', sample_col=None, group_col=None, label_row=None, transpose=True) → AnnData`. Pass the actual separator your file uses (`'__'` for the malaria demo, `'/'` for many MS-DIAL exports).

Peak annotation:
- `ov.metabol.annotate_peaks(mz, *, polarity='positive'|'negative', ppm=10.0, custom_adducts=None, mass_db=None) → pd.DataFrame` with `peak_idx`, `mz`, `adduct`, `kegg`, `name`, `mass_diff_ppm`. Returns one row per (peak, candidate compound).
- `ov.metabol.mummichog_basic(mz, pvalue, *, polarity='positive', ppm=10.0, significance_cutoff=0.05, n_perm=1000, min_overlap=2, pathways=None, mass_db=None, seed=0) → pd.DataFrame` with `pathway`, `overlap`, `set_size`, `pvalue`, `padj`. Pure-Python port — input is the `mz` array and the per-peak `pvalue` array (not a DEG DataFrame).
- `ov.metabol.mummichog_external(mz, pvalue, retention_time=None, *, mode='pos'|'neg', significance_cutoff=0.05, outdir=None, **kwargs) → pd.DataFrame` — thin wrapper over Li's reference `mummichog` PyPI package; install with `pip install mummichog` for full LIBSDB / activity-network features. **Polarity uses `mode='pos'/'neg'` here** (not `polarity='positive'/'negative'` like `mummichog_basic` accepts) because the upstream `Mummichog` class names it that way. `**kwargs` is forwarded into the upstream package's config, so extras the `mummichog` CLI accepts (e.g. `instrument`, `ppm`, network-scoring knobs) can be passed through there.
- `ov.metabol.load_pathways() → dict[str, list[str]]` — `{pathway_name: [kegg_id, ...]}`. Default fetches full KEGG (~550 pathways), cached under `~/.cache/omicverse/metabol/`.
- `ov.metabol.fetch_chebi_compounds() → pd.DataFrame` with `chebi`, `name`, `formula`, `mw` (monoisotopic), `kegg`, `hmdb`. Cached after first call (~54k compounds). Pass back into `annotate_peaks(mass_db=...)` / `mummichog_basic(mass_db=...)` to avoid re-fetching.
- `ov.metabol.fetch_kegg_pathways(organism=None) → dict` — direct KEGG REST fetch.

Lipidomics ingest + parsing:
- `ov.metabol.parse_lipid(name: str) → Optional[LipidIdentity]` — returns `None` for unparseable names so the caller can filter.
- `LipidIdentity` dataclass: `lipid_class`, `total_carbons`, `total_db`, `backbone`, `raw`; convenience methods `is_saturated()` and `is_polyunsaturated(threshold=2)`.
- `ov.metabol.annotate_lipids(adata, *, feature_names=None) → AnnData` (returns a copy; adds `var['lipid_class']`, `var['total_carbons']`, `var['total_db']`).
- `ov.metabol.aggregate_by_class(adata, *, agg='sum'|'mean') → AnnData` — class-level totals matrix.
- `ov.metabol.lion_enrichment(hits, background, *, ontology=None, min_size=3) → pd.DataFrame` — over-representation across LION lipid-ontology terms.
- `ov.metabol.fetch_lion_associations() → dict` — full LION ontology (cached); pass to `lion_enrichment(ontology=...)`.

Recognized lipid classes (canonical, extensible — see `LIPID_CLASSES` in `_lipidomics.py`):
- Glycerophospholipids: `PC, PE, PS, PG, PI, PA` and lyso variants `LPC, LPE, LPS, LPG, LPI, LPA`.
- Sphingolipids: `SM, Cer, GlcCer, LacCer, Hex2Cer, Hex3Cer`.
- Glycerolipids: `TAG, TG, DAG, DG, MAG, MG`.
- Sterols / fatty acids: `CE, FA, Chol, BMP`.

## Boundary

**Inside scope:**
- LC-MS ingest with m/z + RT parsing (`read_lcms`).
- m/z + p-value → mummichog (basic and external backends).
- Per-peak adduct annotation against ChEBI/KEGG.
- Lipid name parsing (LIPID MAPS shorthand sum-composition level).
- Lipid-class aggregation matrix.
- LION ontology enrichment.

**Outside scope — separate skill:**
- Imputation / normalization / transformation (preprocessing skill — but **don't impute LC-MS zeros**, that's enforced upstream).
- HMDB / KEGG / ChEBI ID mapping for *named* metabolites (`ov.metabol.map_ids` — pathway-multifactor skill).
- MSEA ORA / GSEA on named metabolites with KEGG / Reactome / WikiPathways (pathway-multifactor skill).
- Multivariate PLS-DA / OPLS-DA / VIP — the multivariate skill works fine on lipidomics matrices once `aggregate_by_class` (optional) and Pareto-scale (mandatory) have run.
- Targeted lipidomics with full sn-1/sn-2 species resolution — current parser is sum-composition only.

## Branch Selection

**`mummichog_basic` vs. `mummichog_external`**
- `mummichog_basic` — pure-Python port; ships with KEGG-only pathway DB; reproducible with `seed`. Use for portable runs and CI.
- `mummichog_external` — wraps the upstream `mummichog` PyPI package (Li 2013 reference implementation). Adds full adduct table, activity-network scoring, LIBSDB / SMPDB pathways. Use when you need publication-grade multi-DB enrichment and have `pip install mummichog` available.

**`use_pvalue=True` on `volcano` for untargeted LC-MS**: with thousands of peaks vs. <20 samples, BH-FDR effectively zeros out everything. The honest visualization is raw p-value; mummichog's permutation null then absorbs the FDR burden at the *pathway* level. Set `clip_log2fc=8.0` (or similar) so a few below-detection zeros don't blow the x-axis.

**`ppm` for adduct matching**: `10.0` is the typical Orbitrap / TOF tolerance. Tighten to 5 ppm for high-resolution Orbitrap; loosen to 20 ppm for older Q-TOF data. Wider ppm = more candidates per peak but also more false positives.

**`polarity`** must match the ionization mode of the run. Mixed-polarity datasets need split processing (positive ions through `polarity='positive'`, negative through `'negative'`, then merge results) — there's no auto-detect.

**`significance_cutoff` for mummichog**: default 0.05 picks peaks at raw p<0.05 as the "hit" set. For very small cohorts (n<10/group), consider 0.10. The *pathway* p-values are still calibrated by `n_perm`, so loosening the per-peak cutoff doesn't double-count.

**`n_perm` for mummichog**: 1000 for reportable runs (gives 0.001-resolution pathway p-values). For interactive exploration 200–500 is fine. The cost is linear: each permutation re-shuffles peak labels and re-tabulates pathway hits.

**`min_overlap` for mummichog / LION**: 2 is the floor (single-overlap pathways are noise); 3 is canonical for publication. Higher means fewer but more reliable hits.

**Lipidomics — class vs. species level**
- Run `differential` on the species-level (cleaned) matrix when you need fine-grained chain-length / saturation differences.
- Run on the `aggregate_by_class` matrix when you want class-level summaries (e.g. "PC up, PE down") — a parsimonious model with fewer multiple-testing corrections.
- LION enrichment expects species-level *names* in the hit list, not aggregated classes — feed it `deg.index` from the species-level test.

**Caching** — for any LC-MS run that calls mummichog more than once (e.g. positive + negative), explicitly fetch `pathways = ov.metabol.load_pathways()` and `mass_db = ov.metabol.fetch_chebi_compounds()` once, then pass `pathways=pathways, mass_db=mass_db` into each call. The functions auto-fetch on first use but the network round-trip is slow.

## Input Contract

- LC-MS ingest: file is a wide CSV with feature names that encode m/z and RT separated by `feature_id_sep` (e.g. `212.0345__3.45` for `212.0345 m/z, 3.45 RT`). The reader writes `var['m_z']` and `var['rt']` as floats; `obs[group]` from `label_row` if provided.
- For `mummichog_basic`/`_external`: `mz` and `pvalue` are 1-D arrays of equal length. `mz` is in Daltons, not Th. Bad m/z values (NaN / non-positive) are silently dropped — confirm `len(mz) == adata.n_vars` before passing.
- Lipidomics ingest: `var_names` must look like `'PC 34:1'`, `'Cer d18:1/24:0'`, `'TAG 54:3'`, `'LPC 18:0'`, `'TAG 54:3;O'` (the trailing `;O` for oxidized species is recognized). Anything outside the `LIPID_CLASSES` registry parses to `NaN` and is dropped.
- After `aggregate_by_class`, the class-level AnnData carries the original `obs` only if you copy it manually (`class_adata.obs = adata.obs`); the function focuses on `X` and `var` only.

## Minimal Execution Patterns

```python
# Untargeted LC-MS — full path
import omicverse as ov
import numpy as np
ov.plot_set()

adata = ov.metabol.read_lcms(
    'malaria_feature_table.csv',
    feature_id_sep='__',
    label_row='Label',
    transpose=True,
)

# Preprocess: PQN + log (DON'T impute LC-MS zeros)
adata_p = ov.metabol.normalize(adata, method='pqn')
adata_p = ov.metabol.transform(adata_p, method='log')

# Differential — per-peak p-values are the mummichog input
deg = ov.metabol.differential(
    adata_p, group_col='group',
    group_a='Semi_immue', group_b='Naive',
    method='welch_t', log_transformed=True,
)

# Volcano on raw p-values + clipped log2fc — the honest untargeted view
fig, ax = ov.metabol.volcano(deg, padj_thresh=0.01, log2fc_thresh=2.0,
                             label_top_n=0, use_pvalue=True, clip_log2fc=8.0)

# Annotate peaks against KEGG via [M+H]+ adducts
ann = ov.metabol.annotate_peaks(adata.var['m_z'].values,
                                polarity='positive', ppm=10.0)

# Mummichog pathway enrichment — pre-fetch DBs to avoid re-network
pathways = ov.metabol.load_pathways()
mass_db = ov.metabol.fetch_chebi_compounds()
mumm = ov.metabol.mummichog_basic(
    mz=adata.var['m_z'].values,
    pvalue=deg['pvalue'].values,
    polarity='positive', ppm=10.0,
    significance_cutoff=0.05,
    n_perm=1000, min_overlap=2,
    mass_db=mass_db, pathways=pathways,
    seed=0,
)
mumm.head(10)[['pathway', 'overlap', 'set_size', 'pvalue', 'padj']]
```

```python
# Lipidomics — full path
import pandas as pd, anndata as ad

# Build AnnData (replace with your loader; brca demo uses raw CSVs)
matrix = pd.read_csv('brca_matrix.csv', index_col=0)
clin   = pd.read_csv('brca_clin.csv',  index_col=0)
adata = ad.AnnData(matrix.T)             # samples × lipids
adata.obs = clin.loc[adata.obs_names].copy()
adata.obs['group'] = adata.obs['Group']  # rename to canonical 'group'

# Sanity-check the parser
for name in ['PC 34:1', 'TAG 54:3', 'Cer d18:1/24:0', 'LPC 18:0']:
    print(name, '→', ov.metabol.parse_lipid(name))

# Annotate the matrix
adata = ov.metabol.annotate_lipids(adata)
print(adata.var['lipid_class'].value_counts())

# Class-level totals (optional)
class_adata = ov.metabol.aggregate_by_class(adata, agg='sum')
class_adata.obs = adata.obs

# Differential at the species level (drop unparseable first)
adata_clean = adata[:, adata.var['lipid_class'].notna()].copy()
deg = ov.metabol.differential(
    adata_clean, group_col='group',
    group_a='Cancer', group_b='Benign',
    method='welch_t', log_transformed=True,
)
fig, ax = ov.metabol.volcano(deg, padj_thresh=0.10, log2fc_thresh=0.3, label_top_n=8)

# LION enrichment from species-level hits
hits = deg[deg.padj < 0.10].index.tolist()
background = deg.index.tolist()
lion = ov.metabol.lion_enrichment(hits, background, min_size=2)
lion.head(10)[['term', 'category', 'overlap', 'set_size',
               'odds_ratio', 'pvalue', 'padj']]
```

## Validation

- After `read_lcms`: `'m_z' in adata.var.columns` and `'rt' in adata.var.columns`; both are floats. If not, `feature_id_sep` was wrong.
- Before mummichog: `len(adata.var['m_z']) == len(deg)` and the orderings match (the function consumes positionally — feed `adata.var['m_z'].values` and `deg.loc[adata.var_names, 'pvalue'].values` if `deg` may have been re-ordered).
- Mummichog returning an empty frame is a *result*, not an error — log it with the cohort size and `min_overlap` so the user knows to relax thresholds.
- Sanity-test mummichog by feeding a synthesized hit set drawn from a known KEGG pathway (TCA in the tutorial); the seeded pathway should appear at the top.
- After `annotate_lipids`: `var['lipid_class'].notna().sum()` should account for most of the matrix; if half the matrix is `NaN`, the data isn't using LIPID MAPS shorthand and you need a different ingest.
- After `aggregate_by_class`: `class_adata.n_vars` equals the number of distinct `lipid_class` values; `class_adata.X.sum(axis=1)` should equal `adata_clean.X.sum(axis=1)` when `agg='sum'` (within float tolerance).
- LION enrichment: report empirical p-values with `min_size>=2` only; single-member hits are uninformative.

## Resource Map

- See [`reference.md`](reference.md) for copy-paste-ready snippets per workflow.
- See [`references/source-grounding.md`](references/source-grounding.md) for verified `_mummichog` / `_lipidomics` / `_fetchers` signatures and the `LipidIdentity` docstring backfill log.
- For ID mapping of *named* metabolites and MSEA pathway enrichment, see the `omicverse-bulk-metabol-pathway-multifactor` skill.
- For preprocessing chain and LC-MS-specific drift / SERRF / ComBat, see the `omicverse-bulk-metabol-preprocessing` skill.
- For multivariate / biomarker analysis on lipidomics class- or species-level matrices, see the `omicverse-bulk-metabol-multivariate` skill.

## Examples
- "Load `malaria_feature_table.csv` LC-MS data with `feature_id_sep='__'`, run PQN+log, do differential for Semi_immue vs Naive, and run mummichog on the m/z + pvalue with 1000 permutations."
- "Annotate peaks against KEGG via [M+H]+ adducts at 10 ppm, then sanity-check mummichog by seeding TCA-cycle compounds with low p-values and confirming TCA shows up at the top."
- "Parse the `brca_matrix.csv` lipidomics matrix, annotate every lipid species, drop unparseable, run Cancer-vs-Benign differential, then LION-enrich the padj<0.10 hits."
- "Aggregate the lipidomics matrix to class-level totals and propagate the clinical metadata."

## References
- Tutorial notebooks:
  - [`t_metabol_04_untargeted.ipynb`](https://omicverse.readthedocs.io/en/latest/Tutorials-metabol/t_metabol_04_untargeted/) — malaria plasma LC-MS, mummichog.
  - [`t_metabol_05_lipidomics.ipynb`](https://omicverse.readthedocs.io/en/latest/Tutorials-metabol/t_metabol_05_lipidomics/) — `lipidr` BRCA demo, LION.
- Live API verified against `omicverse.metabol` — see [`references/source-grounding.md`](references/source-grounding.md).
