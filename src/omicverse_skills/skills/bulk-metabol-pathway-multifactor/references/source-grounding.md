# Source Grounding — Pathway, Multifactor, DGCA, MOFA, MTBLS1

## Interfaces Checked

`omicverse.metabol` v0.4.0 (master `78586cbe`+) plus `omicverse.utils.load_metabolights`. Verified via `inspect.signature` + `inspect.getdoc` + direct reading of `omicverse/metabol/{_msea,_multifactor,_correlation,_integration,_id_mapping,plotting}.py` and the metabolights ingest helper. Cell-by-cell cross-checking against `t_metabol_03_pathway`, `t_metabol_07_multifactor`, `t_metabol_09_dgca`, `t_metabol_10_multiomics`, `t_metabol_11_real_data_mtbls1`.

## Live signatures — MSEA & ID mapping

```python
ov.metabol.msea_ora(
    hits: Iterable[str],
    background: Iterable[str], *,
    pathways: Optional[dict[str, list[str]]] = None,
    min_size: int = 3,
    mass_db: Optional[pd.DataFrame] = None,
) -> pd.DataFrame  # pathway, overlap, set_size, odds_ratio, pvalue, padj

ov.metabol.msea_gsea(
    deg: pd.DataFrame, *,
    stat_col: str = 'stat',
    pathways: Optional[dict[str, list[str]]] = None,
    n_perm: int = 1000,
    min_size: int = 3,
    max_size: int = 500,
    seed: int = 0,
    mass_db: Optional[pd.DataFrame] = None,
) -> pd.DataFrame  # es, nes, pval, fdr, matched_size, genes, ledge_genes (pathway in INDEX)

ov.metabol.map_ids(
    names: Iterable[str], *,
    targets: tuple[str, ...] = ('hmdb', 'kegg', 'chebi'),
    mass_db: pd.DataFrame | None = None,
) -> pd.DataFrame
```

Source-grounded behavior:
- `msea_ora` uses `scipy.stats.fisher_exact` per pathway; multiple-testing via BH-FDR. Internally name → KEGG ID via `map_ids`.
- `msea_gsea` wraps `gseapy.prerank` (vendored under `omicverse.external.gseapy`). Output column names follow vendored gseapy lowercase convention; pathway term is the DataFrame **index**.
- `stat_col='stat'` is the canonical signed test statistic from `differential`. Setting `'log2fc'` ranks by effect size; `'pvalue'` ranks by raw p (less standard but supported).
- `map_ids(targets=...)` accepts any subset of `('hmdb', 'kegg', 'chebi', 'pubchem', 'lipidmaps')`.
- Resolution coverage on the shipped lookup is partial — passing `mass_db=ov.metabol.fetch_chebi_compounds()` widens it via formula + monoisotopic-mass fuzzy matching.

## Live signatures — multi-factor designs

```python
ov.metabol.asca(
    adata, *,
    factors: Sequence[str],
    include_interactions: bool = True,
    n_components: int = 2,
    n_permutations: int = 0,
    layer: Optional[str] = None,
    center: bool = True,
    seed: int = 0,
) -> ASCAResult

ov.metabol.mixed_model(
    adata, *,
    formula: str,
    groups: str,
    re_formula: Optional[str] = '1',
    term: Optional[str] = None,
    layer: Optional[str] = None,
) -> pd.DataFrame  # per-feature: coef, stderr, tvalue, pvalue, padj, ci_low, ci_high

ov.metabol.meba(
    adata, *,
    group_col: str,
    time_col: str,
    subject_col: str,
    groups: Optional[tuple] = None,
    layer: Optional[str] = None,
) -> pd.DataFrame  # per-feature: F, pvalue, padj
```

Source-grounded behavior:
- `asca` partitions total SS by ANOVA on each factor, runs SVD on each effect's mean matrix. `n_permutations=500` is canonical for reportable runs; `0` skips significance testing.
- `ASCAResult` exposes `.effects: dict[str, ASCAEffect]` keyed by factor name (interactions named `'A:B'`); `.summary()` is one-row-per-effect; `.residual_ss`, `.total_ss` for variance partition; per-effect `.scores`, `.loadings`, `.variance_explained`.
- `mixed_model` runs `statsmodels.MixedLM.from_formula` per feature; the `term` argument selects which fixed-effect coefficient to report (e.g. `'treatment[T.drug]'` for a categorical contrast — patsy `T.<level>` syntax).
- `re_formula='1'` is the random-intercept default; `'1 + time'` would add a random slope. Forwards to `MixedLM.from_formula(re_formula=...)`.
- `meba` implements the Tai-Speed 2006 Hotelling T² test on per-subject time-trajectories — detects features whose temporal *pattern* differs between groups (interaction-style signal). The tutorial validates on synthetic data where features 0..4 carry interaction signal and 5..9 carry only main effects; MEBA recovers the planted interaction features.

## Live signatures — differential correlation

```python
ov.metabol.dgca(
    adata, *,
    group_col: str,
    group_a: Optional[str] = None,
    group_b: Optional[str] = None,
    features: Optional[Sequence[str]] = None,
    method: str = 'pearson',
    abs_r_threshold: float = 0.3,
    layer: Optional[str] = None,
) -> pd.DataFrame
# columns: feature_a, feature_b, r_a, r_b, p_a, p_b,
#          z_diff, p_diff, padj, dc_class

ov.metabol.corr_network(
    adata, *,
    group_col: Optional[str] = None,
    group: Optional[str] = None,
    method: str = 'pearson',
    abs_r_threshold: float = 0.3,
    padj_threshold: Optional[float] = None,
    layer: Optional[str] = None,
) -> pd.DataFrame
# columns: feature_a, feature_b, r, pvalue, padj
# .attrs['n_samples'] records the cohort size
```

Source-grounded behavior:
- `dgca`'s `z_diff` is the Fisher-z-transformed difference between the two correlations; `p_diff` is the standard normal p; `padj` is BH-FDR across all retained pairs.
- `dc_class` is a string encoding the rewiring direction: `'+/+'`, `'-/-'`, `'+/-'`, `'-/+'`, `'+/0'`, `'0/+'`, `'-/0'`, `'0/-'`, `'0/0'`. The order is `class_a / class_b`; `0` means correlation didn't pass `abs_r_threshold` in that group.
- `abs_r_threshold` filters at the *pair* level (retain if `max(|r_a|, |r_b|) >= threshold`), so reversal pairs `+/-` and `-/+` are still retained when one side is barely above threshold.
- `corr_network` retains only pairs passing both `abs_r_threshold` and (optionally) `padj_threshold`; the per-condition cohort size is stored in `attrs['n_samples']` so the user can spot underpowered networks.

## Live signatures — multi-omics MOFA

```python
ov.metabol.run_mofa(
    views: dict, *,
    n_factors: int = 10,
    outfile: str | Path = 'mofa_model.hdf5',
    scale_views: bool = True,
    center_groups: bool = True,
    max_iter: int = 500,
    convergence_mode: str = 'fast',
    gpu_mode: bool = False,
    seed: int = 0,
) -> pd.DataFrame  # samples × retained_factors
```

Source-grounded behavior:
- Bridges to `mofapy2` (lazy-imported). Each view is a sample-aligned `AnnData`; the function asserts `obs_names` match across views in the same order.
- `scale_views=True` z-scores each view independently — necessary when views have different units.
- The trained MOFA model persists to `outfile` (HDF5); the returned DataFrame is the post-pruning factor matrix — MOFA+ drops factors whose variance contribution falls below an internal threshold (typically 0.01).
- `convergence_mode` is forwarded to `mofapy2`: `'fast'` for prototyping, `'medium'` / `'slow'` for reportable runs.

## Live signatures — Metabolights ingest

```python
ov.utils.load_metabolights(
    study_id: str, *,
    group_col: Optional[str] = None,
    cache_dir: str | Path = 'metabolights_cache',
    maf_name: Optional[str] = None,
    sample_name_col: str = 'Sample Name',
    refresh: bool = False,
) -> AnnData
```

Source-grounded behavior:
- Downloads the ISA-Tab files + MAF from `ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/<study_id>/` (HTTPS-mirrored). Caches under `cache_dir`; `refresh=True` bypasses cache.
- Joins the MAF (one row per metabolite, columns include `metabolite_identification`, `m/z`, `retention time`, `database_identifier`) with the sample sheet (one row per sample, factor columns `Factor Value[<name>]`).
- When `group_col` is given, the column is renamed to `'group'` in `obs` to match the rest of `ov.metabol`'s convention. Common `Factor Value` columns: `'Metabolic syndrome'`, `'Gender'`, `'Age'`, etc. — the exact form is `'Factor Value[<name>]'` because that's how ISA-Tab stores them.
- The MTBLS1 study contains `metabolite_identification` as a `var` column; the case study uses it for human-readable VIP labels and pathway-name mapping.
- Multi-MAF studies (e.g. positive- + negative-mode runs): pass `maf_name='m_LCMS_pos.tsv'` to disambiguate.

## Docstring supplementation log

| Symbol | Prior state | Action |
|---|---|---|
| `metabol.plotting.corr_network_plot` | 1-line | filled — full Parameters block (layout, node_size, edge_width_scale, r_column, with_labels, label_font_size, seed); empty-edges placeholder behavior; semantics of edge colour (red = positive r, blue = negative). |
| `metabol.plotting.asca_variance_bar` | 1-line | filled — semantics of the residual bar as model-fit sanity check; description of variance-fraction percentage labels. |

`msea_ora` (29L), `msea_gsea` (26L), `map_ids` (23L), `asca` (37L), `mixed_model` (32L), `meba` (34L), `dgca` (41L), `corr_network` (31L), `run_mofa` (38L), `pathway_bar` (20L), `pathway_dot` (17L), `dgca_class_bar` (already expanded in preprocessing skill), `load_metabolights` (full Numpy-style) — all already well-documented; left as-is. `ASCAResult` had a 4-line class docstring with usage example; left as-is (the per-effect attributes are typed via `ASCAEffect`).

## Notebook ↔ skill alignment

| Notebook | Mapped to |
|---|---|
| `t_metabol_03` `differential` + `map_ids` + `msea_ora` + `pathway_bar` / `pathway_dot` + `msea_gsea` (column-name pitfall) + VIP × DEG cross-check | Quick Workflow §1-5 |
| `t_metabol_07` `asca` with permutation null + `mixed_model` with patsy formula + `meba` time-series + synthetic-recovery validation | Quick Workflow §6-8; Validation block |
| `t_metabol_09` `dgca` + `dgca_class_bar` + per-condition `corr_network` + set arithmetic on edge keys | Quick Workflow §9-12 |
| `t_metabol_10` `run_mofa` with sample-aligned dual views | Quick Workflow §13-15 |
| `t_metabol_11` `load_metabolights` + sample QC + cv_filter + impute / normalize / log / pareto + differential + volcano + opls_da + s_plot + vip_bar (with `metabolite_identification`) + msea_ora + asca over phenotype × sex + roc_feature + biomarker_panel + dgca + corr_network + corr_network_plot | Quick Workflow §16-18; full case-study block in `reference.md` |

## Reviewer-Run Empirical Checks

- All cited functions importable: `from omicverse.metabol import msea_ora, msea_gsea, map_ids, asca, mixed_model, meba, dgca, corr_network, corr_network_plot, dgca_class_bar, asca_variance_bar, run_mofa, ASCAResult, pathway_bar, pathway_dot` ✓; `from omicverse.utils import load_metabolights` ✓.
- Confirmed `msea_gsea` returns DataFrame with pathway in **index** (vendored gseapy convention); skill explicitly calls out the `reset_index()` requirement before `pathway_dot`.
- Confirmed `dgca` `dc_class` uses the 9-state `+/-/0` × `+/-/0` encoding; `dgca_class_bar` palette matches.
- Confirmed `corr_network.attrs['n_samples']` carries the per-condition cohort size (not the parent `adata.n_obs`).
- No live smoke run executed; the synthetic-recovery validation in `t_metabol_07` (MEBA) and `t_metabol_03`-style ID-coverage check are documented as in-skill acceptance probes.
