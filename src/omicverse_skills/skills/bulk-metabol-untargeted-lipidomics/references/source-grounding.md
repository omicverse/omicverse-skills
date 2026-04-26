# Source Grounding — Untargeted LC-MS & Lipidomics

## Interfaces Checked

`omicverse.metabol` v0.4.0 (`master` `78586cbe`+). Verified via `inspect.signature` + `inspect.getdoc` + direct reading of `omicverse/metabol/{_lipidomics,_mummichog,_fetchers,io}.py`. Cell-by-cell cross-checking against `t_metabol_04_untargeted.ipynb` and `t_metabol_05_lipidomics.ipynb`.

## Live signatures — LC-MS ingest

```python
ov.metabol.read_lcms(
    path: str | Path, *,
    feature_id_sep: str = '/',
    sample_col: Optional[str] = None,
    group_col: Optional[str] = None,
    label_row: Optional[str] = None,
    transpose: bool = True,
) -> AnnData
```

Source-grounded behavior:
- Splits each var name on `feature_id_sep`; first token → `var['m_z']` (float, Da), second token → `var['rt']` (float, min).
- **Note the column name is `m_z` (with underscore), not `mz`** — earlier draft of the preprocessing skill mis-spelled this; downstream code must use `adata.var['m_z'].values`.
- `label_row` selects a row of the file as the group factor (the malaria demo uses `'Label'`); written into `obs['group']`.
- `transpose=True` (the default for LCMS exports) flips rows-as-features files into the AnnData `samples × features` convention.
- Zeros in the matrix are kept as zeros (LC-MS convention for below-detection); they are *not* converted to NaN.

## Live signatures — peak annotation & mummichog

```python
ov.metabol.annotate_peaks(
    mz: np.ndarray, *,
    polarity: str = 'positive',
    ppm: float = 10.0,
    custom_adducts: Optional[list[tuple[str, float, str]]] = None,
    mass_db: Optional[pd.DataFrame] = None,
) -> pd.DataFrame  # columns: peak_idx, mz, adduct, kegg, name, mass_diff_ppm

ov.metabol.mummichog_basic(
    mz: np.ndarray,
    pvalue: np.ndarray, *,
    polarity: str = 'positive',
    ppm: float = 10.0,
    significance_cutoff: float = 0.05,
    n_perm: int = 1000,
    min_overlap: int = 2,
    pathways: Optional[dict[str, list[str]]] = None,
    mass_db: Optional[pd.DataFrame] = None,
    seed: int = 0,
) -> pd.DataFrame  # pathway, overlap, set_size, pvalue, padj

ov.metabol.mummichog_external(
    mz: np.ndarray,
    pvalue: np.ndarray, *,
    polarity: str = 'positive',
    ppm: float = 10.0,
    significance_cutoff: float = 0.05,
    ...
) -> pd.DataFrame
```

Source-grounded behavior:
- `mummichog_basic` consumes the **raw `pvalue` array, not `padj`** — the canonical mummichog formulation expects unadjusted per-peak p-values. The function then runs its own permutation null at the *pathway* level, which is what makes it usable on small-n untargeted runs where BH-FDR zeros everything.
- Peaks with `pvalue < significance_cutoff` form the hit set; the rest are background. The two sets are passed through adduct-resolved KEGG mapping using the same `mass_db` as `annotate_peaks`.
- `n_perm` permutations shuffle the peak labels (hit/non-hit) and rebuild the pathway hit table; the permutation null is per-pathway. p-value resolution is `1/n_perm`.
- `min_overlap` filters pathways with too few mapped hits — pathways with `overlap=1` are never reported.
- Both `pathways` and `mass_db` are optional kwargs; passing them avoids re-fetching from KEGG / ChEBI on every call. The DBs cache to `~/.cache/omicverse/metabol/`.
- `mummichog_external` is a thin wrapper over `mummichog` PyPI package; raises `ImportError` if not installed.

## Live signatures — fetchers / pathway DB

```python
ov.metabol.load_pathways(...) -> dict[str, list[str]]
ov.metabol.fetch_kegg_pathways(organism: Optional[str] = None) -> dict
ov.metabol.fetch_chebi_compounds(...) -> pd.DataFrame  # chebi, name, formula, mw, kegg, hmdb
ov.metabol.fetch_lion_associations(...) -> dict
ov.metabol.fetch_hmdb_from_name(...) -> dict
ov.metabol.clear_cache(...) -> None
```

Source-grounded behavior:
- `load_pathways` returns the full KEGG pathway database by default (~550 pathways). First call hits `https://rest.kegg.jp`; subsequent calls read from disk cache.
- `fetch_chebi_compounds` downloads ChEBI flat-file TSVs (`compounds.tsv.gz`, `chemical_data.tsv.gz`, `database_accession.tsv.gz`) over HTTPS from EBI's FTP, joins them, and returns a master compound table (~54k rows). The `mw` column is monoisotopic mass; `kegg` and `hmdb` columns are joined cross-references.
- `fetch_lion_associations` returns `{term_name: {"category": str, "members": [lipid_class, ...]}}` — same shape as the bundled `lion_subset.json` so `lion_enrichment` consumes either format.
- All fetchers cache responses; `clear_cache()` resets if the upstream sources updated.

## Live signatures — lipidomics

```python
ov.metabol.parse_lipid(name: str) -> Optional[LipidIdentity]
ov.metabol.annotate_lipids(adata, *, feature_names=None) -> AnnData
ov.metabol.aggregate_by_class(adata, *, agg: str = 'sum') -> AnnData
ov.metabol.lion_enrichment(
    hits: Iterable[str],
    background: Iterable[str], *,
    ontology: Optional[dict[str, dict]] = None,
    min_size: int = 3,
) -> pd.DataFrame

@dataclass
class LipidIdentity:
    lipid_class: str
    total_carbons: int
    total_db: int
    backbone: Optional[str] = None
    raw: str = ""
    def is_saturated() -> bool: ...
    def is_polyunsaturated(threshold: int = 2) -> bool: ...
```

Source-grounded behavior:
- The `LIPID_CLASSES` registry in `_lipidomics.py` lists 27 canonical classes; the regex matches longest-first so `LPC` doesn't accidentally match the `PC` class. To extend: add to `LIPID_CLASSES` + recompile `_PATTERN`.
- `parse_lipid` returns `None` for unparseable names — caller filters/annotates as "not a lipid". `annotate_lipids` propagates this as `var['lipid_class'] = NaN`.
- `annotate_lipids` returns a *copy* of `adata` (existing var columns preserved); adds `lipid_class`, `total_carbons`, `total_db` columns.
- `aggregate_by_class` requires `var['lipid_class']` to exist; result has `n_vars = n_lipid_classes` (NaN class is dropped). `agg='sum'` is canonical for relative-quantification lipidomics; `agg='mean'` is sometimes useful for absolute-quantification data.
- `lion_enrichment` uses Fisher's exact + BH-FDR; `min_size` filters single-member terms. Categories include lipid class, fatty-acid composition, function, and physical properties.

## Docstring supplementation log

| Symbol | Prior state | Action |
|---|---|---|
| `LipidIdentity` (class) | 1-line | filled — added Attributes block (lipid_class / total_carbons / total_db / backbone / raw with semantics), explanation of sum-composition level, convenience methods documented. |
| `LipidIdentity.is_saturated` | empty | filled (1 line — "True if no double bonds — typical SFA-rich species"). |
| `LipidIdentity.is_polyunsaturated` | empty | filled (1 line — "True if degree of unsaturation ≥ threshold (default 2 — PUFA)"). |

`parse_lipid` (4 lines), `annotate_lipids` (5 lines), `aggregate_by_class` (6 lines) had short but accurate docstrings — left as-is; the class-level `LipidIdentity` docstring now carries the deeper semantic context that those functions reference.

## Notebook ↔ skill alignment

| Notebook section | Skill section |
|---|---|
| `t_metabol_04` ingest with `feature_id_sep='__'`, `label_row='Label'` | Quick Workflow §1 (untargeted) |
| `t_metabol_04` PQN + log preprocessing | Demoted — handled by preprocessing skill (with the LC-MS-zero rule called out) |
| `t_metabol_04` differential + raw-pvalue volcano with `clip_log2fc` | Quick Workflow §3-4 |
| `t_metabol_04` `annotate_peaks` adduct lookup | Quick Workflow §6 |
| `t_metabol_04` `mummichog_basic` with `n_perm=1000`, pre-fetched DBs | Quick Workflow §5; Branch Selection (basic vs external) |
| `t_metabol_04` synthetic-pathway sanity check (TCA seeding) | Quick Workflow §7; reference.md sanity-check block |
| `t_metabol_05` LIPID MAPS sanity-check (`parse_lipid` on demo strings) | Quick Workflow §2 (lipidomics) |
| `t_metabol_05` `annotate_lipids` + class distribution | Quick Workflow §3 |
| `t_metabol_05` `aggregate_by_class` with metadata propagation | Quick Workflow §4 |
| `t_metabol_05` species-level differential + volcano | Quick Workflow §5 |
| `t_metabol_05` `lion_enrichment` on hits | Quick Workflow §6 |

## Reviewer-Run Empirical Checks

- All cited functions importable: `from omicverse.metabol import read_lcms, annotate_peaks, mummichog_basic, mummichog_external, load_pathways, fetch_chebi_compounds, fetch_lion_associations, parse_lipid, annotate_lipids, aggregate_by_class, lion_enrichment, LipidIdentity` ✓
- Confirmed `var['m_z']` (with underscore) is the column name; corrected the lookup in skill text.
- Tutorial cells in both notebooks use exactly the verified signatures.
- No live smoke run executed (per project policy); all evidence is signature- and source-level. The `mummichog_basic` synthetic-TCA recovery test in the tutorial cell IS the closest thing to a smoke check and is documented in `reference.md`.
