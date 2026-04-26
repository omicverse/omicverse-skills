# Source Grounding — 16S amplicon end-to-end

## Interfaces Checked

`omicverse.alignment` (vsearch + dada2 backends + amplicon AnnData composer) and `omicverse.micro` v0.x (master `78586cbe`+). Verified via `inspect.signature` + `inspect.getdoc` + direct reading of `omicverse/alignment/{amplicon_16s_pipeline,vsearch,dada2,build_amplicon_anndata}.py` and `omicverse/micro/{_diversity,_ord,_da,_pp,_phylo}.py`. Cell-by-cell cross-checking against `t_16s_amplicon.ipynb` and `t_16s_dada2.ipynb`.

## Live signatures — pipeline + ingest

```python
ov.alignment.amplicon_16s_pipeline(
    fastq_dir: Optional[str] = None,
    samples: Optional[Sequence[Tuple[str, str, Optional[str]]]] = None,
    workdir: Optional[str] = None,
    db_fasta: Optional[str] = None, *,
    primer_fwd: Optional[str] = None, primer_rev: Optional[str] = None,
    backend: str = 'vsearch',          # 'vsearch' | 'dada2'
    threads: int = 4,
    jobs: Optional[int] = None,
    merge_max_diffs: int = 10, merge_min_overlap: int = 16,
    filter_max_ee: float = 1.0, filter_min_len: int = 0, filter_max_len: int = 0,
    derep_min_uniq: int = 2,
    unoise_alpha: float = 2.0, unoise_minsize: int = 2,
    chimera_removal: bool = True,
    otutab_identity: float = 0.97,
    sintax_cutoff: float = 0.8, sintax_strand: str = 'both',
    sample_metadata: Optional[pd.DataFrame] = None,
    overwrite: bool = False,
) -> AnnData

ov.alignment.fetch_rdp(db_dir: Optional[str] = None, overwrite: bool = False) -> str

ov.alignment.build_amplicon_anndata(
    otutab_tsv: str,
    asv_fasta: str,
    sintax_tsv: Optional[str] = None,
    sample_metadata: Optional[pd.DataFrame] = None,
    sample_order: Optional[Sequence[str]] = None,
) -> AnnData
```

Source-grounded behavior:
- `amplicon_16s_pipeline` is a fluent wrapper: when `backend='vsearch'`, dispatches to `vsearch.merge_pairs / filter_quality / dereplicate / unoise3 / uchime3_denovo / sintax / usearch_global` then calls `build_amplicon_anndata`. When `backend='dada2'`, dispatches to the `dada2` submodule (same logical chain, different implementation).
- `fastq_dir` and `samples` are mutually exclusive — pass one or the other.
- `workdir` is mandatory and never falls back to `$HOME` or `/tmp` (deliberate choice; documented in the `workdir` parameter doc).
- `primer_fwd` / `primer_rev`: if both set, runs `cutadapt`; otherwise skips primer trimming. The mothur SOP test dataset has primers pre-removed, so the tutorial uses `None`.
- `chimera_removal=True` runs `uchime3_denovo` after UNOISE3; canonical for 16S amplicon.
- `sample_metadata`: `pd.DataFrame` indexed by sample-id; written into `adata.obs`. Sample-id mismatch ⇒ silently dropped from final AnnData (caller should validate counts match).
- The returned `AnnData` has a sparse integer count matrix in `.X`, ASV centroid sequences as a `var` column, and 7-rank SINTAX taxonomy in `var` (`domain, phylum, class, order, family, genus, species, sintax_confidence`). Empty strings for ranks below the bootstrap cutoff.

## Live signatures — `ov.micro` diversity / ordination / DA

```python
class Alpha:
    def __init__(self, adata, rarefy_depth: Optional[int] = None, seed: int = 0): ...
    def run(self, metrics=('shannon', 'observed_otus'), write_to_obs: bool = True,
            tree_key: str = 'tree') -> pd.DataFrame: ...
    def shannon(self) -> pd.Series: ...
    def observed(self) -> pd.Series: ...

class Beta:
    def __init__(self, adata, rarefy_depth: Optional[int] = None, seed: int = 0): ...
    def run(self, metric: str = 'braycurtis', rarefy: Optional[bool] = None,
            tree_key: str = 'tree', write_to_obsp: bool = True) -> pd.DataFrame: ...
    def braycurtis(self, rarefy: bool = True) -> pd.DataFrame: ...

class Ordinate:
    def __init__(self, adata, dist_key: str = 'braycurtis'): ...
    def pcoa(self, n: int = 3, write_to_obsm: bool = True) -> pd.DataFrame: ...
    def nmds(self, n: int = 2, random_state: int = 0,
             write_to_obsm: bool = True) -> pd.DataFrame: ...
    def proportion_explained(self) -> Optional[np.ndarray]: ...

class DA:
    def __init__(self, adata): ...
    def wilcoxon(self, group_key, group_a=None, group_b=None,
                 rank=None, relative=True, min_prevalence=0.1) -> pd.DataFrame: ...
    def deseq2(self, group_key, group_a=None, group_b=None,
               rank=None, min_prevalence=0.1, alpha=0.05) -> pd.DataFrame: ...
    def ancombc(self, group_key, rank=None,
                min_prevalence=0.1, pseudocount=1.0) -> pd.DataFrame: ...

# Module-level helpers
ov.micro.rarefy(adata, depth=None, seed=0, drop_shallow=True,
                save_original=True, copy=False) -> AnnData
ov.micro.filter_by_prevalence(adata, min_prevalence=0.1, min_count=1, copy=False) -> AnnData
ov.micro.collapse_taxa(adata, rank='genus', unassigned_label='Unassigned') -> AnnData
ov.micro.clr(adata, layer_out='clr', copy=False) -> AnnData
ov.micro.ilr(adata, layer_out='ilr', copy=False) -> AnnData
```

Source-grounded behavior:
- `Alpha` rarefies to `rarefy_depth` if set, otherwise to the shallowest sample. Internally uses `skbio.diversity.alpha_diversity`. Faith PD requires a tree at `adata.uns[tree_key]`.
- `Beta.run`'s `rarefy=None` defaults to `True` for `braycurtis`/`jaccard`/`unifrac_*`, `False` for `aitchison`. UniFrac requires `unifrac` package; CLR / Aitchison runs without scikit-bio.
- `Ordinate.pcoa` uses `skbio.stats.ordination.pcoa`, which returns eigenvalues + variance-explained; the variance-explained vector is stored at `adata.uns['micro'][f'{dist_key}_pcoa_var']` and accessed via `proportion_explained()`.
- `Ordinate.nmds` uses `sklearn.manifold.MDS(dissimilarity='precomputed', normalized_stress='auto')`; final stress at `adata.uns['micro'][f'{dist_key}_nmds_stress']`.
- `DA.wilcoxon` runs Mann-Whitney per feature; `relative=True` divides by the per-sample sum (so it's testing on relative abundances, not raw counts). `rank=None` keeps ASV-level; otherwise calls `collapse_taxa` first.
- `DA.deseq2` requires `pip install pydeseq2`; raises `ImportError` otherwise. Runs on raw counts directly (no relative-scaling pre-step).
- `DA.ancombc` requires `skbio>=0.7.1`; pseudocount added to zero-counts before composition transform.
- `rarefy(save_original=True)` caches `adata.X` at `adata.layers['raw_counts']` so subsequent transforms don't lose the original counts.

## Docstring supplementation log

Driven by this skill's source-grounding pass:

| Symbol | Prior state | Action |
|---|---|---|
| `Alpha.shannon` | empty | filled — semantics + depth-rarefaction warning |
| `Alpha.observed` | empty | filled — depth-dependent caveat |
| `Beta.braycurtis` | empty | filled — Bray-Curtis interpretation, rarefy default |
| `Ordinate.nmds` | 1-line | filled — Parameters block; rank-preservation rationale; obsm/uns persistence |
| `Ordinate.proportion_explained` | 1-line | filled — None-when-NMDS caveat; downstream usage |
| `DA.deseq2` | 1-line | filled — full Parameters block; method positioning vs Wilcoxon / ANCOM-BC; output schema |

`MMvec.fit / cooccurrence / conditional_probabilities / top_pairs / plot_mmvec_training / plot_embedding_biplot` were also expanded as part of the same module audit; those land in the micro-metabol-paired skill's source-grounding doc rather than this one (different skill scope), but the diff is committed in `omicverse/micro/_pair.py`.

Other public symbols (`Alpha.run`, `Beta.run`, `Ordinate.pcoa`, `DA.wilcoxon`, `DA.ancombc`, `rarefy`, `filter_by_prevalence`, `collapse_taxa`, `clr`, `ilr`, `attach_tree`, `combine_studies`, `meta_da`, `paired_spearman`, `paired_cca`, `simulate_paired`, `fetch_franzosa_ibd_2019`) — all already had ≥10-line Numpy-style docstrings; left as-is.

`omicverse.alignment.dada2.{merge_pairs, make_seqtab, remove_chimeras}` are flagged by the registry as missing docstrings, but these are stepwise dada2 internals not directly exposed in the tutorial — the public path is `amplicon_16s_pipeline(backend='dada2')` which has a complete docstring. Out of scope for this skill's docstring backfill; tracked for a future `alignment.dada2` audit.

## Notebook ↔ skill alignment

| Notebook section | Skill section |
|---|---|
| `t_16s_amplicon` `fetch_rdp` + `amplicon_16s_pipeline` (one-shot) | Quick Workflow §1-3, Minimal Execution Patterns |
| `t_16s_amplicon` stepwise vsearch chain (merge → filter → derep → unoise3 → uchime3 → sintax → usearch_global → build_amplicon_anndata) | Interface Summary; Minimal Execution Patterns (stepwise block) |
| `t_16s_amplicon` Alpha (`shannon`, `observed_otus`, `simpson`) at min depth | Quick Workflow §5; Validation |
| `t_16s_amplicon` Beta (Bray-Curtis) + Ordinate.pcoa + group-coloured PC1/PC2 scatter | Quick Workflow §6; reference.md PCoA block |
| `t_16s_amplicon` `DA.wilcoxon` at phylum, Early vs Late | Quick Workflow §7 |
| `t_16s_dada2` same FASTQs with `backend='dada2'`, `filter_max_ee=2.0`; recovers same diversity / taxonomy patterns | Branch Selection (vsearch vs dada2) |
| `t_16s_dada2` phylum stacked-bar + Shannon/Observed boxplots + DADA2 PCoA | reference.md (delegates to plotting helpers) |

## Reviewer-Run Empirical Checks

- All cited functions importable: `from omicverse.alignment import amplicon_16s_pipeline, fetch_rdp, build_amplicon_anndata; from omicverse.alignment import vsearch; from omicverse.micro import Alpha, Beta, Ordinate, DA, rarefy, filter_by_prevalence, collapse_taxa, clr, ilr, attach_tree` ✓
- `ov.micro.Alpha(adata).run(metrics=['shannon', 'observed_otus'])` writes both columns into `adata.obs` (verified by reading the source).
- `ov.micro.Ordinate(adata, dist_key='braycurtis').pcoa(n=3)` writes `adata.obsm['braycurtis_pcoa']` and `adata.uns['micro']['braycurtis_pcoa_var']` (verified by reading the source — both notebooks read these exact keys).
- No live smoke run executed; the mothur SOP demo dataset is the established acceptance vehicle for both tutorials and would serve as the future smoke target.
