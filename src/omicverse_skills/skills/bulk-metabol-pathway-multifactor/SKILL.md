---
name: omicverse-bulk-metabol-pathway-multifactor
description: Pathway interpretation, multi-factor designs, differential correlation, MOFA multi-omics, and the MTBLS1 real-data case study on a preprocessed metabolomics AnnData. Use when you need MSEA ORA / GSEA pathway enrichment with KEGG IDs, multi-factor ANOVA-SCA / mixed-model / MEBA designs, DGCA differential correlation networks, MOFA+ joint factorization across metabolomics + RNA-seq, or a worked end-to-end T2D study from MetaboLights ingestion to biomarker panel.
---

# OmicVerse Bulk Metabolomics — Pathway, Multifactor, DGCA & MOFA

## Goal

Cover the four downstream interpretation paths that consume a preprocessed metabolomics `AnnData` and produce a biological story:

1. **Pathway / MSEA** — resolve metabolite names to external IDs (HMDB / KEGG / ChEBI / PubChem / LIPID MAPS) and run over-representation (ORA) or GSEA-style ranked enrichment, plot bars / dots.
2. **Multi-factor designs** — ASCA (ANOVA-SCA), mixed models (MixedLM), MEBA time-series — for studies with treatment × time × patient layouts.
3. **Differential correlation (DGCA)** — find metabolite *pairs* whose correlation rewires between conditions, plus per-condition correlation networks.
4. **Multi-omics MOFA+** — joint factorization across metabolomics and RNA-seq (or any other AnnData-aligned views).

Plus a worked end-to-end real-data case study: **MTBLS1 (urine NMR, Type 2 Diabetes)** that ties preprocessing → univariate / multivariate / pathway / multifactor / DGCA / network into one runbook.

## Quick Workflow

### Pathway / MSEA (`t_metabol_03`)

1. Run `differential` on the preprocessed (PQN+log) `AnnData` to get a per-metabolite p-value DataFrame.
2. (Optional) Resolve metabolite names to external IDs with `ov.metabol.map_ids(...)` to confirm KEGG coverage.
3. Run `ov.metabol.msea_ora(hits, background, min_size=3)` for over-representation, or `ov.metabol.msea_gsea(deg, stat_col='stat', n_perm=500)` for GSEA-style ranked enrichment.
4. Plot the result: `pathway_bar(...)` for raw p-value bars, `pathway_dot(...)` for the canonical dot plot (size = overlap count, x = odds_ratio or NES, color = p-value).
5. Cross-check by fitting OPLS-DA, pulling the VIP table, and intersecting `padj<0.20 ∧ VIP>1` for a tighter shortlist.

### Multi-factor (`t_metabol_07`)

6. **ASCA** — for crossed factor designs (e.g. treatment × time): `ov.metabol.asca(adata, factors=['treatment', 'time'], include_interactions=True, n_permutations=500)`. Returns an `ASCAResult` with per-effect variance and permutation p-values; `result.summary()` is a one-row-per-effect DataFrame; `asca_variance_bar(result)` plots the variance partition.
7. **Mixed models** — for repeated measures with random subject effect: `ov.metabol.mixed_model(adata, formula='treatment + time', groups='patient', term='treatment[T.drug]')`. Per-feature MixedLM fit; uses patsy formula syntax; the `term` argument selects which fixed-effect coefficient to report.
8. **MEBA** — two-group time-series (Hotelling T² across timepoints, Tai 2006): `ov.metabol.meba(adata, group_col='group', time_col='time', subject_col='subject')`. Detects features whose temporal *trajectory* differs between groups (interaction-style signal).

### DGCA / correlation networks (`t_metabol_09`)

9. Run `ov.metabol.dgca(adata, group_col='group', group_a, group_b, method='spearman'|'pearson', abs_r_threshold=0.3)` for differential correlation between conditions. Returns pairs with `r_a`, `r_b`, `z_diff`, `padj`, `dc_class` (the rewiring class — `+/+`, `-/-`, `+/-`, `-/+`, `+/0`, ...).
10. Plot the class breakdown with `ov.metabol.dgca_class_bar(dc)` — log-scaled counts; reversals (`+/-`, `-/+`) are red.
11. For per-condition networks, run `ov.metabol.corr_network(adata, group_col='group', group=<label>, method='spearman', abs_r_threshold=0.5, padj_threshold=0.05)` once per condition, then compare edge sets (shared / condition-specific).
12. Visualize with `ov.metabol.corr_network_plot(edges_df, layout='spring', label_font_size=6)`. For dense graphs, slice the top-N edges first.

### Multi-omics MOFA (`t_metabol_10`)

13. Build sample-aligned AnnData per view (metabol view + RNA view, **same `obs_names` in same order**). Each view is its own `AnnData`.
14. Run `factors = ov.metabol.run_mofa(views={'metabol': adata_m, 'rna': adata_r}, n_factors=5, outfile='mofa.hdf5', scale_views=True, max_iter=200, seed=0)`. Returns the factor matrix (samples × retained factors) — MOFA+ prunes inactive factors automatically.
15. The trained model is persisted to `outfile` (HDF5); load it with the standard MOFA+ API to get per-view loadings.

### MTBLS1 case study (`t_metabol_11`)

16. Ingest with `ov.utils.load_metabolights('MTBLS1', group_col='Factor Value[Metabolic syndrome]', cache_dir='mtbls1_demo')` — auto-downloads from EBI's Metabolights public mirror, joins the MAF with the sample sheet, and writes `obs['group']`.
17. Sample-QC (`sample_qc` + `sample_qc_plot`) → CV-filter (`cv_filter(across='all', cv_threshold=1.5)`) → impute → normalize → log → Pareto.
18. Run the full downstream: differential + volcano, PLS-DA + OPLS-DA + S-plot + VIP, MSEA ORA, ASCA over `group + Gender`, ROC + biomarker panel, DGCA + correlation network — using human-readable metabolite names from `var['metabolite_identification']` for plotting.

## Interface Summary

ID mapping (already discussed in untargeted-lipidomics skill for m/z-based annotation; here for *name*-based ID resolution):
- `ov.metabol.map_ids(names, *, targets=('hmdb', 'kegg', 'chebi'), mass_db=None) → pd.DataFrame` indexed by name with columns for each requested target.
- `ov.metabol.normalize_name(name) → str` — strip parenthetical synonyms, lowercase, dedupe spaces (used internally; expose for diagnostic queries).

MSEA pathway enrichment:
- `ov.metabol.msea_ora(hits, background, *, pathways=None, min_size=3, mass_db=None) → pd.DataFrame` — Fisher's exact test against KEGG (default) or any user-supplied `{pathway: [kegg_id, ...]}` dict.
- `ov.metabol.msea_gsea(deg, *, stat_col='stat', pathways=None, n_perm=1000, min_size=3, max_size=500, seed=0, mass_db=None) → pd.DataFrame` — GSEA-style ranked enrichment via `gseapy.prerank`. Column names follow vendored gseapy: `es`, `nes`, `pval`, `fdr`, `matched_size`, `genes`, `ledge_genes`. Pathway term is the **DataFrame index, not a column** — `reset_index()` before passing to `pathway_dot`.
- `ov.metabol.load_pathways() → dict[str, list[str]]` (full KEGG pathway database, ~550 pathways, cached on first call).

Multi-factor designs:
- `ov.metabol.asca(adata, *, factors, include_interactions=True, n_components=2, n_permutations=0, layer=None, center=True, seed=0) → ASCAResult`. ASCA = ANOVA on each effect → SVD on each effect's mean matrix → variance fraction + per-feature loadings.
- `ASCAResult.effects: dict[str, ASCAEffect]` — per-effect data; `.summary() → pd.DataFrame` (one row per effect: `variance_explained`, `permutation_pvalue`, `n_components`); `.residual_ss`, `.total_ss`.
- `ov.metabol.mixed_model(adata, *, formula, groups, re_formula='1', term=None, layer=None) → pd.DataFrame` — per-feature `statsmodels.MixedLM` fit. `formula` is patsy syntax (e.g. `'treatment + time'`); `groups` is the random-intercept column; `term` selects the fixed-effect coefficient to report (use `coefname[T.level]` for categorical contrasts).
- `ov.metabol.meba(adata, *, group_col, time_col, subject_col, groups=None, layer=None) → pd.DataFrame` — per-feature Hotelling T² for two-group time-courses (Tai 2006). Detects features whose temporal pattern differs between groups.

DGCA + correlation networks:
- `ov.metabol.dgca(adata, *, group_col, group_a=None, group_b=None, features=None, method='pearson'|'spearman', abs_r_threshold=0.3, layer=None) → pd.DataFrame` — pair-level differential correlation; `dc_class` encodes rewiring direction.
- `ov.metabol.corr_network(adata, *, group_col=None, group=None, method='pearson'|'spearman', abs_r_threshold=0.3, padj_threshold=None, layer=None) → pd.DataFrame` — within-condition pairwise correlations passing thresholds; `attrs['n_samples']` records the cohort size used.
- `ov.metabol.corr_network_plot(edges_df, *, layout='spring'|'circular'|'kamada_kawai', node_size=70, edge_width_scale=2.5, r_column='r', label_font_size=7, ...)`.
- `ov.metabol.dgca_class_bar(dc_df, *, log=True, ax=None, figsize=(6, 3.5))` — log-scaled bar chart of DC-class counts.

Plotting helpers (pathway):
- `ov.metabol.pathway_bar(enrichment, *, term_col=None, score_col='pvalue', top_n=15, ax=None, figsize=(6, 5))`.
- `ov.metabol.pathway_dot(enrichment, *, term_col=None, size_col='overlap', x_col='odds_ratio', color_col='pvalue', top_n=15, ax=None, figsize=(6.5, 5.5))`.
- `ov.metabol.asca_variance_bar(asca_result, *, ax=None, figsize=(5.5, 3))`.

Multi-omics MOFA+:
- `ov.metabol.run_mofa(views: dict, *, n_factors=10, outfile='mofa_model.hdf5', scale_views=True, center_groups=True, max_iter=500, convergence_mode='fast', gpu_mode=False, seed=0) → pd.DataFrame`. Returns the post-pruning factor matrix; the full model lives in `outfile`.

MetaboLights ingest:
- `ov.utils.load_metabolights(study_id, *, group_col=None, cache_dir='metabolights_cache', maf_name=None, sample_name_col='Sample Name', refresh=False) → AnnData`.

## Boundary

**Inside scope:**
- All of MSEA (ORA + GSEA) on named metabolites, with KEGG (default) or any user-supplied pathway DB.
- ASCA, mixed-model, MEBA — the three multi-factor design helpers.
- DGCA, `corr_network`, network plotting.
- MOFA+ joint factorization with metabolomics + any aligned omics view.
- Metabolights ingest + the MTBLS1 worked example.

**Outside scope — separate skill:**
- m/z-based mummichog enrichment for *untargeted* LC-MS — see the untargeted-lipidomics skill.
- LION ontology enrichment for lipid species — same skill above.
- Imputation / normalization / transformation / drift / SERRF / ComBat / sample-outlier QC — preprocessing skill.
- PLS-DA / OPLS-DA fitting itself (the case study uses them; the skill *referenced* is `omicverse-bulk-metabol-multivariate`, not this one).
- Biomarker panel + nested CV — multivariate skill.
- Drug-perturbation Cox / signature scoring — out of scope for `ov.metabol`.

## Branch Selection

**MSEA — ORA vs GSEA**
- ORA when you have a clear hit/non-hit cutoff (e.g. `padj<0.20`) and a moderate background (~50–500 metabolites). Outputs an odds ratio + p-value per pathway.
- GSEA when the signal is distributed across many small effects rather than concentrated in a few big hits. Pass the full DEG table with `stat_col='stat'` (or `'log2fc'`) — no cutoff. Outputs NES + permutation p / FDR per pathway. The vendored `gseapy.prerank` puts pathway names in the **index**; use `reset_index()` before plotting.
- Sample sizes too small for either? Skip MSEA — small-cohort metabolomics has too few hits for either test to calibrate.

**`map_ids` target selection**
- `('hmdb', 'kegg', 'chebi')` is the canonical trio; `'pubchem'` and `'lipidmaps'` are also supported.
- Coverage on the shipped lookup is partial — for the cachexia 63-metabolite NMR set, most resolve; for novel LC-MS features, expect 60–80 % unresolved. Pre-fetch `mass_db=ov.metabol.fetch_chebi_compounds()` and pass it to widen the resolver.

**Multi-factor — ASCA vs MixedLM vs MEBA**
- `asca` for **crossed factor designs with replicates** (e.g. treatment × time × subject), where you want a global variance partition plus per-feature loadings. Requires balanced data; permutation p-values calibrate effect significance.
- `mixed_model` for **repeated measures** (same subjects measured in multiple conditions) — provides per-feature p-values for any fixed-effect coefficient. Random subject intercept handles within-subject correlation.
- `meba` for **two-group time-courses** where the question is "which features behave differently *over time* between groups". Single-test, Hotelling T²-based.
- Three groups with no time? Use `anova` from the preprocessing skill.

**DGCA — Spearman vs Pearson**
- Spearman is the default for metabolomics: it tolerates the non-linear, non-Gaussian distributions typical of NMR / MS data.
- Pearson only when the features are confirmed approximately log-normal (after log transform, this is reasonable for many platforms — but verify with QQ plots before trusting Pearson p-values).
- `abs_r_threshold` filters weak correlations on *both* sides — pairs are retained only if `max(|r_a|, |r_b|) >= threshold`. Tighter threshold = fewer but stronger rewiring candidates.

**`corr_network` thresholds**
- For network *visualization*: `abs_r_threshold=0.5, padj_threshold=0.05` keeps the figure readable.
- For network *enumeration / shared-edge analysis*: `abs_r_threshold=0.3, padj_threshold=0.05` retains more structure for set arithmetic.
- Always pair `abs_r_threshold` with `padj_threshold` — a high `|r|` on small n can be noise.

**MOFA — `n_factors` and `scale_views`**
- `n_factors=10` is a reasonable starting point — MOFA+ prunes inactive factors during fitting, so over-asking is cheap.
- `scale_views=True` (default) z-scores each view independently — necessary when views have different units (intensity vs counts vs scaled).
- `convergence_mode='fast'` for prototyping; `'medium'` / `'slow'` for reportable runs.
- `gpu_mode=True` only if you've installed `mofapy2` with CUDA support — typically not needed for the small cohorts metabolomics studies use.

**Metabolights ingest**
- `group_col='Factor Value[<phenotype>]'` — the canonical Metabolights factor column convention.
- `cache_dir` writes to disk; subsequent runs reuse cached files. Set `refresh=True` to bypass.
- `maf_name` only needed when a study ships multiple MAFs (e.g. positive- and negative-mode runs).

## Input Contract

- For all four workflows, start from a preprocessed `AnnData`: `obs['group']` populated, `X` is **log-transformed**. Use a Pareto-scaled copy for ASCA (which is variance-partitioning and benefits from Pareto), but keep a log-only copy for `differential` / `roc_feature` / `dgca` / `mixed_model`.
- For MOFA: every view's `AnnData` has identical `obs_names` in the same order. The function will raise on misalignment.
- For Metabolights ingest: an internet connection (first call) and write access to `cache_dir`. Subsequent calls work offline.
- For ASCA / `mixed_model`: every factor in `factors`/`formula` must exist as an `obs` column with the right dtype (categorical for ASCA factors; categorical or numeric for MixedLM).
- For `corr_network`: avoid running on the full 500+ metabolite matrix without thresholding — the pair count grows quadratically.

## Minimal Execution Patterns

```python
# MSEA — ORA + GSEA on the cachexia data
import omicverse as ov
import matplotlib.pyplot as plt
ov.plot_set()

adata = ov.metabol.read_metaboanalyst('human_cachexia.csv', group_col='Muscle loss')
adata = ov.metabol.normalize(adata, method='pqn')
adata = ov.metabol.transform(adata, method='log')

deg = ov.metabol.differential(adata, group_col='group',
                              group_a='cachexic', group_b='control',
                              method='welch_t', log_transformed=True)

# Optional ID coverage check
ids = ov.metabol.map_ids(adata.var_names.tolist())
print(f'{(ids["kegg"] != "").sum()} / {len(ids)} metabolites resolve to KEGG')

# ORA
hits = deg[deg.padj < 0.20].index.tolist()
background = deg.index.tolist()
ora = ov.metabol.msea_ora(hits, background, min_size=3)
ov.metabol.pathway_bar(ora, score_col='pvalue', top_n=10)
ov.metabol.pathway_dot(ora, size_col='overlap', x_col='odds_ratio',
                       color_col='pvalue', top_n=10)

# GSEA — pathway names are in the INDEX
gsea = ov.metabol.msea_gsea(deg, stat_col='stat', n_perm=500, seed=0)
gsea_plot = gsea.reset_index().rename(columns={gsea.index.name or 'index': 'pathway'})
ov.metabol.pathway_dot(gsea_plot, size_col='matched_size', x_col='nes',
                       color_col='pval', top_n=10)
plt.show()
```

```python
# Multi-factor — ASCA, mixed model, MEBA
asca_res = ov.metabol.asca(
    adata, factors=['treatment', 'time'],
    include_interactions=True,
    n_components=2, n_permutations=500, seed=0,
)
print(asca_res.summary())
ov.metabol.asca_variance_bar(asca_res); plt.show()

mlm = ov.metabol.mixed_model(
    adata,
    formula='treatment + time',
    groups='patient',
    term='treatment[T.drug]',
)
print(mlm.sort_values('pvalue').head())

# Time-series: needs subject_col for repeated measures within group
meba_tbl = ov.metabol.meba(
    adata_ts, group_col='group',
    time_col='time', subject_col='subject',
)
print(meba_tbl.sort_values('pvalue').head())
```

```python
# DGCA + per-condition correlation networks
dc = ov.metabol.dgca(adata, group_col='group',
                     group_a='cachexic', group_b='control',
                     method='spearman', abs_r_threshold=0.3)
ov.metabol.dgca_class_bar(dc); plt.show()

edges_a = ov.metabol.corr_network(adata, group_col='group', group='cachexic',
                                  method='spearman', abs_r_threshold=0.5,
                                  padj_threshold=0.05)
edges_b = ov.metabol.corr_network(adata, group_col='group', group='control',
                                  method='spearman', abs_r_threshold=0.5,
                                  padj_threshold=0.05)

a_keys = set(edges_a.apply(lambda r: frozenset((r['feature_a'], r['feature_b'])), axis=1))
b_keys = set(edges_b.apply(lambda r: frozenset((r['feature_a'], r['feature_b'])), axis=1))
print(f'shared edges: {len(a_keys & b_keys)}')
print(f'A-only:       {len(a_keys - b_keys)}')

ov.metabol.corr_network_plot(edges_a.head(40), figsize=(7, 6), label_font_size=6)
plt.show()
```

```python
# Multi-omics MOFA — both views must share obs_names in the same order
factors = ov.metabol.run_mofa(
    views={'metabol': adata_m, 'rna': adata_r},
    n_factors=5,
    outfile='mofa_demo.hdf5',
    scale_views=True,
    max_iter=200,
    seed=0,
)
print(f'retained {factors.shape[1]} factors after MOFA+ pruning')
print(factors.head())
```

```python
# MTBLS1 case study — full pipeline
adata = ov.utils.load_metabolights(
    'MTBLS1',
    group_col='Factor Value[Metabolic syndrome]',
    cache_dir='mtbls1_demo',
)

qc = ov.metabol.sample_qc(adata, n_components=3, alpha=0.95)
ov.metabol.sample_qc_plot(qc); plt.show()

adata = ov.metabol.cv_filter(adata, across='all', cv_threshold=1.5)
adata = ov.metabol.impute(adata, method='qrilc', seed=0)
adata = ov.metabol.normalize(adata, method='pqn')
adata_log = ov.metabol.transform(adata, method='log')
adata_pareto = ov.metabol.transform(adata_log, method='pareto', stash_raw=False)

# Univariate + volcano on log-only
deg = ov.metabol.differential(adata_log, group_col='group',
                              group_a='diabetes mellitus', group_b='Control Group',
                              method='welch_t', log_transformed=True)
ov.metabol.volcano(deg, padj_thresh=0.05, log2fc_thresh=0.5, label_top_n=8)

# Multivariate on Pareto
opls = ov.metabol.opls_da(adata_pareto, group_col='group',
                          group_a='diabetes mellitus', group_b='Control Group',
                          n_ortho=1)
mb_names = adata_pareto.var['metabolite_identification'].values
ov.metabol.vip_bar(opls, mb_names, top_n=15)

# MSEA on names
mass_db = ov.metabol.fetch_chebi_compounds()
hit_names = adata.var.loc[deg[deg['padj']<0.05].index, 'metabolite_identification'].tolist()
bg_names = adata.var['metabolite_identification'].dropna().tolist()
enr = ov.metabol.msea_ora(hit_names, bg_names, mass_db=mass_db)
ov.metabol.pathway_bar(enr, top_n=10)

# Multi-factor: phenotype + sex
asca_res = ov.metabol.asca(adata_pareto,
                           factors=['group', 'Factor Value[Gender]'],
                           include_interactions=True, n_permutations=100, seed=0)
ov.metabol.asca_variance_bar(asca_res)

# DGCA + network
dc = ov.metabol.dgca(adata_log, group_col='group',
                     group_a='diabetes mellitus', group_b='Control Group',
                     method='spearman')
ov.metabol.dgca_class_bar(dc)

edges_t2d = ov.metabol.corr_network(adata_log, group_col='group',
                                    group='diabetes mellitus',
                                    method='spearman',
                                    abs_r_threshold=0.5, padj_threshold=0.05)
named = adata.var['metabolite_identification'].to_dict()
edges_plot = edges_t2d.copy()
edges_plot['feature_a'] = edges_plot['feature_a'].map(named).fillna(edges_plot['feature_a'])
edges_plot['feature_b'] = edges_plot['feature_b'].map(named).fillna(edges_plot['feature_b'])
ov.metabol.corr_network_plot(edges_plot.head(40), figsize=(7, 6), label_font_size=6)
```

## Validation

- After `msea_ora` / `msea_gsea`: an empty result is a *signal*, not an error — log the cohort size and `min_size`/`n_perm` so the user knows whether to relax thresholds. For ORA, also report background size; with a small background (<50 metabolites), the test is underpowered.
- After `msea_gsea`: pathway names are in the **index**, not a column. Code that does `gsea['pathway']` will KeyError; use `gsea.reset_index()`.
- After `asca`: `result.summary()` reports per-effect `variance_explained` and `permutation_pvalue`. A residual variance fraction > 0.5 suggests the model misses dominant signal; revisit the factors. Permutation p < 0.05 is the canonical effect-significance bar.
- After `mixed_model`: the `term` argument must match a coefficient *exactly* — for categorical, that's `coefname[T.level]` (e.g. `'treatment[T.drug]'`); typos return all-NaN results.
- After `meba`: validate on a synthetic interaction-only signal that the planted features rank above non-interaction features (the tutorial does this — replicate before trusting on novel data).
- After `dgca`: `dc_class='+/-'` and `'-/+'` are the strongest rewiring; `'+/+'` and `'-/-'` are concordant (no rewiring); `'0/0'` is no signal in either group. Use the `dgca_class_bar` to gut-check the breakdown.
- After `corr_network`: `attrs['n_samples']` records the per-condition cohort size — use this for power calibration, not `adata.n_obs`.
- After `run_mofa`: confirm `factors.shape[1] < n_factors` (MOFA+ pruned at least one factor) — if not, the model didn't converge on a parsimonious solution; rerun with `convergence_mode='medium'`.
- For the MTBLS1 case study: the log + Pareto matrices serve different downstream stages — passing the wrong one inflates fold-changes (Pareto into `differential`) or degrades `Q²` (log-only into `opls_da`). Keep both copies.

## Resource Map

- See [`reference.md`](reference.md) for compact copy-paste snippets per workflow.
- See [`references/source-grounding.md`](references/source-grounding.md) for verified `_msea` / `_multifactor` / `_correlation` / `_integration` / `_id_mapping` signatures and the `corr_network_plot` / `asca_variance_bar` docstring backfill log.
- For preprocessing chain that produces the input `AnnData`, see `omicverse-bulk-metabol-preprocessing`.
- For PLS-DA / OPLS-DA / VIP / S-plot / biomarker panel cited in the MTBLS1 case study, see `omicverse-bulk-metabol-multivariate`.
- For untargeted m/z-based pathway inference (mummichog) and lipidomics LION enrichment, see `omicverse-bulk-metabol-untargeted-lipidomics`.

## Examples
- "Run MSEA ORA at padj<0.20 on the cachexia DEG, plot a pathway dot plot, and cross-validate with GSEA on the same DEG ranking."
- "Treatment × time crossed design: run ASCA with 500 permutations and plot the variance partition."
- "Repeated-measures study: per-feature MixedLM with `formula='treatment + time'` and `groups='patient'`, report the `treatment[T.drug]` coefficient."
- "Detect metabolite pairs whose Spearman correlation rewires between cachexic and control (DGCA) and plot the per-condition correlation networks."
- "Joint MOFA+ on metabolomics + RNA-seq with 5 factors, scale_views=True, save the model to mofa_demo.hdf5."
- "Reproduce the MTBLS1 worked example end-to-end: ingest, sample-QC, CV-filter, full preprocessing, multivariate, MSEA, ASCA, biomarker panel, DGCA, correlation network."

## References
- Tutorial notebooks:
  - [`t_metabol_03_pathway.ipynb`](https://omicverse.readthedocs.io/en/latest/Tutorials-metabol/t_metabol_03_pathway/) — MSEA ORA + GSEA on cachexia.
  - [`t_metabol_07_multifactor.ipynb`](https://omicverse.readthedocs.io/en/latest/Tutorials-metabol/t_metabol_07_multifactor/) — ASCA, mixed model, MEBA.
  - [`t_metabol_09_dgca.ipynb`](https://omicverse.readthedocs.io/en/latest/Tutorials-metabol/t_metabol_09_dgca/) — DGCA + correlation networks.
  - [`t_metabol_10_multiomics.ipynb`](https://omicverse.readthedocs.io/en/latest/Tutorials-metabol/t_metabol_10_multiomics/) — MOFA+ on synthetic two-view data.
  - [`t_metabol_11_real_data_mtbls1.ipynb`](https://omicverse.readthedocs.io/en/latest/Tutorials-metabol/t_metabol_11_real_data_mtbls1/) — MTBLS1 worked example.
- Live API verified — see [`references/source-grounding.md`](references/source-grounding.md).
