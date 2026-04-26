---
name: omicverse-microbiome-16s-amplicon-dada2
description: 16S rRNA amplicon analysis from raw FASTQs to a samples × ASVs AnnData with 7-rank SINTAX taxonomy, plus the canonical alpha / beta / ordination / DA stack. Use when running the vsearch-or-DADA2 amplicon pipeline (`ov.alignment.amplicon_16s_pipeline`), when ingesting an existing OTU/ASV count table with `build_amplicon_anndata`, or when computing Shannon / Bray-Curtis / PCoA / Wilcoxon-DA on the resulting AnnData.
---

# OmicVerse Microbiome — 16S amplicon end-to-end

## Goal

Take paired-end Illumina FASTQs from a 16S rRNA gene survey through to an analysis-ready `AnnData` (samples × ASVs, with 7-rank SINTAX taxonomy in `var` and the ASV centroid sequence stored alongside), plus the de-facto downstream stack: alpha-diversity (Shannon / observed OTUs / Simpson) → beta-diversity (Bray-Curtis distance matrix) → PCoA / NMDS ordination → differential abundance (Wilcoxon by default; pyDESeq2 / ANCOM-BC available).

Cover **two backends** behind the same wrapper:
- **vsearch / UNOISE3** — the default; fast (single-threaded C, parallelised via `jobs`); slightly more reads retained per sample on noisy datasets.
- **DADA2** — Bayesian denoising; longer to fit but produces tighter ASV-level resolution. Same wrapper, just `backend='dada2'`.

Plus the **stepwise** API for users who want to inspect or substitute individual stages (merge → filter → dereplicate → UNOISE3 → uchime3 → SINTAX → OTU table → AnnData).

## Quick Workflow

1. **Reference DB**: fetch a SINTAX-formatted 16S reference once with `ov.alignment.fetch_rdp(db_dir=...)` — returns the path to the `.fa.gz` you'll feed all subsequent calls.
2. **Sample metadata**: build a `pd.DataFrame` keyed by sample-id with at least the phenotype/factor column (here called `'group'`). The wrapper writes this into `adata.obs`.
3. **One-shot pipeline**: call `ov.alignment.amplicon_16s_pipeline(fastq_dir=..., workdir=..., db_fasta=..., backend='vsearch'|'dada2', sample_metadata=meta, ...)`. Auto-discovers R1/R2 from FASTQ names, runs the full chain, writes intermediates under `workdir`, returns an `AnnData`.
4. **Validate ingest**: confirm `adata.shape == (n_samples, n_ASVs)`, `'phylum' in adata.var.columns` (and the rest of the 7-rank), `int(adata.X.sum())` is reasonable (>50 % of input reads typically retained).
5. **Alpha diversity**: `ov.micro.Alpha(adata, rarefy_depth=min_depth).run(metrics=['shannon', 'observed_otus', 'simpson'])` — rarefies first, writes per-sample metrics into `adata.obs`.
6. **Beta + ordination**: `ov.micro.Beta(adata, rarefy_depth=min_depth).run(metric='braycurtis')` then `ov.micro.Ordinate(adata, dist_key='braycurtis').pcoa(n=3)`. Distance matrix lives at `adata.obsp['braycurtis']`; PCoA coords at `adata.obsm['braycurtis_pcoa']`; variance-explained at `adata.uns['micro']['braycurtis_pcoa_var']`.
7. **Differential abundance**: `ov.micro.DA(adata).wilcoxon(group_key='group', group_a, group_b, rank='phylum'|'genus'|...)`. Returns DataFrame with effect size + BH-FDR.

## Interface Summary

Reference DB:
- `ov.alignment.fetch_rdp(db_dir=None, overwrite=False) → str` — alias for `fetch_sintax_ref('rdp_16s_v18')`. ~6.8 MB. Returns path to `.fa.gz`.
- `ov.alignment.fetch_sintax_ref(db_name, db_dir=None) → str` — switch DB (`'silva138'`, `'gtdb'`, etc.).

Pipeline wrapper (one-shot path):
```python
ov.alignment.amplicon_16s_pipeline(
    fastq_dir: Optional[str] = None,
    samples: Optional[Sequence[(name, fq1, fq2|None)]] = None,
    workdir: Optional[str] = None,
    db_fasta: Optional[str] = None, *,
    primer_fwd: Optional[str] = None, primer_rev: Optional[str] = None,
    backend: str = 'vsearch',          # or 'dada2'
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
) → AnnData
```

Stepwise vsearch helpers (under `ov.alignment.vsearch.*`):
- `merge_pairs(samples, output_dir, max_diffs=10, min_overlap=16, threads=8, jobs=4)`.
- `filter_quality(merge_res, output_dir, max_ee=1.0, threads=8, jobs=4)`.
- `dereplicate(filt_res, output_dir, min_uniq=2, threads=8)`.
- `unoise3(uniques_fasta, output_dir, alpha=2.0, minsize=2, threads=8)`.
- `uchime3_denovo(asv_fasta, output_dir)`.
- `sintax(asv_fasta, db_fasta, output_dir, cutoff=0.8, strand='both', threads=8)`.
- `usearch_global(combined_fastq, asv_fasta, output_dir, identity=0.97, threads=8)`.

Composer:
- `ov.alignment.build_amplicon_anndata(otutab_tsv, asv_fasta, sintax_tsv=None, sample_metadata=None, sample_order=None) → AnnData` — convert pre-computed OTU table + ASV centroids + SINTAX TSV into the canonical AnnData.

Diversity / ordination / DA (`ov.micro`):
- `ov.micro.Alpha(adata, rarefy_depth=None, seed=0).run(metrics=('shannon', 'observed_otus'), write_to_obs=True, tree_key='tree') → pd.DataFrame`. Convenience methods: `.shannon()`, `.observed()` return individual `pd.Series`. Faith PD via `'faith_pd'` requires a phylogenetic tree at `adata.uns[tree_key]` (see the phylogeny skill).
- `ov.micro.Beta(adata, rarefy_depth=None, seed=0).run(metric='braycurtis'|'jaccard'|'aitchison'|'unifrac_unweighted'|'unifrac_weighted', rarefy=None, tree_key='tree', write_to_obsp=True) → pd.DataFrame`. UniFrac requires `unifrac` package + a tree.
- `ov.micro.Ordinate(adata, dist_key='braycurtis').pcoa(n=3, write_to_obsm=True) → pd.DataFrame` (eigen-decomposition of the distance matrix; `proportion_explained()` gives variance fractions). `.nmds(n=2, random_state=0, write_to_obsm=True)` for non-metric MDS.
- `ov.micro.DA(adata).wilcoxon(group_key, group_a=None, group_b=None, rank=None, relative=True, min_prevalence=0.1) → pd.DataFrame` with `log2fc`, `pvalue`, `padj`. Per-feature Mann-Whitney U on relative abundances.
- `ov.micro.DA(adata).deseq2(group_key, group_a=None, group_b=None, rank=None, min_prevalence=0.1, alpha=0.05)` — pyDESeq2 NB-GLM on raw counts.
- `ov.micro.DA(adata).ancombc(group_key, rank=None, min_prevalence=0.1, pseudocount=1.0)` — ANCOM-BC via skbio≥0.7.1.

Preprocessing helpers:
- `ov.micro.rarefy(adata, depth=None, seed=0, drop_shallow=True, save_original=True, copy=False)` — subsample counts; original counts cached at `adata.layers['raw_counts']` when `save_original=True`.
- `ov.micro.filter_by_prevalence(adata, min_prevalence=0.1, min_count=1, copy=False)` — drop rare features.
- `ov.micro.collapse_taxa(adata, rank='genus', unassigned_label='Unassigned')` — sum ASV counts within taxonomic rank.
- `ov.micro.clr(adata, layer_out='clr', copy=False)` / `ilr(...)` — compositional transforms (post pseudo-count).

## Boundary

**Inside scope:**
- Full vsearch / UNOISE3 pipeline.
- Full DADA2 pipeline (same wrapper, `backend='dada2'`).
- Stepwise vsearch composability + `build_amplicon_anndata` from external outputs.
- 7-rank SINTAX taxonomy assignment with the shipped RDP DB (or any SINTAX-format reference).
- Alpha (Shannon / Simpson / Chao1 / observed OTUs) and beta (Bray-Curtis / Jaccard / Aitchison) diversity.
- PCoA / NMDS ordination on a stored distance matrix.
- Wilcoxon DA at any taxonomic rank.

**Outside scope — separate skill:**
- DA method *comparison* across Wilcoxon / pyDESeq2 / ANCOM-BC — see `omicverse-microbiome-da-comparison`.
- Cross-cohort meta-analysis (`combine_studies`, `meta_da`) — see `omicverse-microbiome-meta-analysis`.
- Phylogenetic tree construction + UniFrac + Faith PD — see `omicverse-microbiome-phylogeny`.
- Paired microbiome + metabolomics integration (MMvec, paired Spearman / CCA) — see `omicverse-micro-metabol-paired`.
- Shotgun metagenomics — `ov.alignment` doesn't ship a metagenomics pipeline; this skill is amplicon-only.

## Branch Selection

**`backend='vsearch'` vs `'dada2'`**
- vsearch (default): faster (~minutes for the mothur SOP demo), ships a single binary, parallelises across samples (`jobs`). Use for production runs and CI.
- DADA2: Bayesian denoising; tighter ASV resolution at the cost of longer runtime. Use when noise floor matters (e.g. low-biomass samples) or when reviewers expect DADA2 specifically.
- Both write the same canonical AnnData layout — switching backends doesn't change downstream code.

**Primer trimming (`primer_fwd` / `primer_rev`)**
- Set both to the actual primer sequences for cutadapt to run first.
- Pass `None` only when you've confirmed primers are already trimmed (e.g. the mothur SOP test dataset).
- Mismatch: cutadapt errors out with a useful message — don't silence it.

**`filter_max_ee` (expected error rate)**
- 1.0 — strict (vsearch default; loses ~30 % of low-quality reads).
- 2.0 — DADA2-tutorial default (more permissive); needed on Q20+ MiSeq runs to retain enough reads.
- Empirical: try 1.0 first, drop to 2.0 if fewer than 50 % of input reads survive the filter.

**`unoise_minsize` and `derep_min_uniq`**
- `unoise_minsize=2` (defaults to 2) — UNOISE3 retains an ASV only if it appears at least this many times across the whole derep table. Lower (1) for low-biomass studies; higher (4–8) to suppress sequencing noise on deep cohorts.
- `derep_min_uniq` is the per-sample equivalent before pooling.

**`sintax_cutoff` (taxonomic confidence)**
- 0.8 (default, RDP recommendation) — bootstrap support threshold; ranks below this become empty strings.
- 0.5 — looser; useful on novel environmental samples where genus-level resolution isn't available in any DB.
- Lower than 0.5 — the published "do-not-go-below" floor; results are not reproducible across DB versions.

**`otutab_identity`**
- 0.97 — classical OTU clustering threshold (≈ species-level for 16S V3-V4); use when you specifically want OTU-style buckets, not ASVs.
- For pure ASV workflows, this stage is just used to build the per-sample count matrix from the original reads — keep at 0.97 for compatibility with downstream tools.

**Alpha metrics**
- `shannon` — diversity (counts both richness and evenness); always report this.
- `observed_otus` — pure richness, depth-sensitive (always rarefy first).
- `simpson` — diversity weighted toward dominant taxa.
- `chao1` — depth-corrected richness estimator; reports unobserved-taxa correction.
- `faith_pd` — phylogenetic diversity; needs a tree (phylogeny skill).

**Beta metrics**
- `braycurtis` — abundance-weighted, the de-facto 16S default.
- `jaccard` — presence/absence; complement to Bray-Curtis when interpretation differs.
- `aitchison` — CLR-Euclidean; correct for compositional data, but harder to interpret intuitively.
- `unifrac_*` — phylogenetically aware; needs a tree (phylogeny skill).

**Ordination — PCoA vs NMDS**
- PCoA: linear, fast, gives variance-explained — use when distances are ≈ Euclidean (Aitchison) or for the canonical "PC1/PC2 with %variance" plot.
- NMDS: rank-based, no eigenvalues; less distorted on Bray-Curtis / Jaccard. Report the final stress (in `adata.uns['micro']`).

**DA method**
- `wilcoxon` — non-parametric, fast, no distributional assumptions. Default for moderate cohorts (n>=10/group).
- `deseq2` — NB-GLM on raw counts; better-calibrated than Wilcoxon, more conservative with small counts. Requires `pip install pydeseq2`.
- `ancombc` — bias-corrected ANCOM; closest to the compositional ground-truth but slowest. Requires `skbio>=0.7.1`. See the DA-comparison skill for trade-offs.

## Input Contract

- Pipeline ingest: paired-end Illumina FASTQs with predictable R1/R2 naming (e.g. `<sample>_R1.fastq.gz`); one sample per pair. The wrapper auto-discovers from `fastq_dir`; pass `samples=[(name, fq1, fq2)]` for explicit control.
- `workdir`: writable directory, `>1 GB` free space recommended for moderate cohorts. The wrapper does **not** fall back to `$HOME` or `/tmp` — pass an explicit path.
- `db_fasta`: SINTAX-formatted 16S FASTA (`fetch_rdp` returns a valid one). Skipping (`db_fasta=None`) means no taxonomy in `var`.
- `sample_metadata`: `pd.DataFrame` indexed by sample-id (matching FASTQ-derived names exactly). The phenotype column is conventionally called `'group'` to match `ov.micro.DA`'s defaults.
- After ingest, `adata.X` is a sparse (or dense) integer count matrix — never log-transform before alpha / beta / DA; those functions handle scaling internally.

## Minimal Execution Patterns

```python
# vsearch one-shot pipeline
import omicverse as ov
import pandas as pd

ov.plot_set()

DB_FASTA = ov.alignment.fetch_rdp(db_dir='./db/rdp')

meta = pd.DataFrame({
    'group': ['Early', 'Early', 'Late', 'Late', 'Mock'],
    'day':   [0, 1, 7, 8, 0],
}, index=['F3D0', 'F3D1', 'F3D7', 'F3D8', 'Mock'])

adata = ov.alignment.amplicon_16s_pipeline(
    fastq_dir='./raw/MiSeq_SOP',
    workdir='./run_vsearch',
    db_fasta=DB_FASTA,
    threads=8, jobs=4,
    primer_fwd=None, primer_rev=None,    # SOP is pre-trimmed
    filter_max_ee=1.0,
    unoise_minsize=2,
    sintax_cutoff=0.8,
    sample_metadata=meta,
)
print(adata)                       # samples × ASVs
print(int(adata.X.sum()), 'reads mapped')
```

```python
# DADA2 backend — same wrapper
adata_dada2 = ov.alignment.amplicon_16s_pipeline(
    samples=[('F3D0', 'raw/F3D0_R1.fastq.gz', 'raw/F3D0_R2.fastq.gz'), ...],
    workdir='./run_dada2',
    db_fasta=DB_FASTA,
    backend='dada2',
    primer_fwd=None, primer_rev=None,
    filter_max_ee=2.0,        # DADA2 typically wants more permissive EE
    sintax_cutoff=0.8,
    threads=4,
    sample_metadata=meta,
)
```

```python
# Stepwise (vsearch) — when you need to inspect / substitute a stage
import numpy as np

merge_res = ov.alignment.vsearch.merge_pairs(samples, output_dir='./step/merged',
                                             max_diffs=10, min_overlap=16,
                                             threads=8, jobs=4)
filt_res  = ov.alignment.vsearch.filter_quality(merge_res, output_dir='./step/filtered',
                                                max_ee=1.0, threads=8, jobs=4)
derep     = ov.alignment.vsearch.dereplicate(filt_res, output_dir='./step/derep',
                                             min_uniq=2, threads=8)
unoise    = ov.alignment.vsearch.unoise3(derep['uniques'], output_dir='./step/asv',
                                         alpha=2.0, minsize=2, threads=8)
nochim    = ov.alignment.vsearch.uchime3_denovo(unoise['asv'], output_dir='./step/asv')
tax       = ov.alignment.vsearch.sintax(nochim['asv'], db_fasta=DB_FASTA,
                                        output_dir='./step/taxonomy',
                                        cutoff=0.8, strand='both', threads=8)
otutab    = ov.alignment.vsearch.usearch_global(derep['combined'], nochim['asv'],
                                                output_dir='./step/otutab',
                                                identity=0.97, threads=8)

adata_step = ov.alignment.build_amplicon_anndata(
    otutab_tsv=otutab['otutab'],
    asv_fasta=nochim['asv'],
    sintax_tsv=tax['tsv'],
    sample_metadata=meta,
    sample_order=[s[0] for s in samples],
)
```

```python
# Diversity + ordination + DA on the AnnData
min_depth = int(np.asarray(adata.X.sum(axis=1)).min())

# alpha — writes shannon/observed_otus/simpson into adata.obs
ov.micro.Alpha(adata, rarefy_depth=min_depth).run(
    metrics=['shannon', 'observed_otus', 'simpson'],
)

# beta — distance matrix in adata.obsp['braycurtis']
ov.micro.Beta(adata, rarefy_depth=min_depth).run(metric='braycurtis')

# PCoA — coords in adata.obsm['braycurtis_pcoa'], var-fraction in uns
ord_ = ov.micro.Ordinate(adata, dist_key='braycurtis').pcoa(n=3)
pct = ord_.proportion_explained() * 100.0

# DA — Wilcoxon at phylum level
da = ov.micro.DA(adata).wilcoxon(
    group_key='group', group_a='Early', group_b='Late',
    rank='phylum', min_prevalence=0.1,
)
print(da.head(10))
```

## Validation

- After ingest: `adata.shape` matches the expected `(n_samples, n_ASVs)`. The mothur SOP demo lands at ~22 samples × ~250 ASVs with vsearch defaults.
- Read retention: `int(adata.X.sum()) / int(<input read total>) > 0.5` is healthy. Below 0.3 indicates the filter / chimera / derep stages dropped too much; revisit `filter_max_ee` and `unoise_minsize`.
- Taxonomy coverage: `(adata.var['phylum'] != '').sum() / adata.n_vars` should exceed 0.8 — lower means the DB is wrong species range or `sintax_cutoff` is too strict.
- Alpha rarefaction: pass `rarefy_depth=int(adata.X.sum(axis=1).min())` to rarefy down to the shallowest sample. If that's <1000, the cohort is too shallow for trustworthy alpha — flag this rather than silently downsampling.
- Beta rarefaction: same — uneven depths inflate Bray-Curtis distances against the deeper samples.
- After ordination: `proportion_explained()` is decreasing; PC1 + PC2 typically capture 30–60 % on 16S. If PC1 alone captures >70 %, suspect a single dominant outlier sample.
- After DA: `(da['padj'] < 0.05).sum()` should be sane for the cohort — for the mothur SOP 22-sample demo, expect ~3–8 phylum-level hits at padj<0.05.
- The DADA2 backend produces *different* ASV identities than vsearch; downstream interpretations (e.g. specific ASV names) don't transfer between backends. Diversity / ordination / DA *patterns* should agree.

## Resource Map

- See [`reference.md`](reference.md) for compact one-shot + stepwise + diversity snippets.
- See [`references/source-grounding.md`](references/source-grounding.md) for verified pipeline + `ov.micro` API signatures and the docstring backfill log (`Alpha.shannon` / `Alpha.observed` / `Beta.braycurtis` / `Ordinate.nmds` / `Ordinate.proportion_explained` / `DA.deseq2`).
- For DA method comparison (Wilcoxon vs pyDESeq2 vs ANCOM-BC), see `omicverse-microbiome-da-comparison`.
- For cross-cohort meta-analysis on the resulting AnnDatas, see `omicverse-microbiome-meta-analysis`.
- For phylogenetic tree + UniFrac + Faith PD, see `omicverse-microbiome-phylogeny`.
- For paired microbiome + metabolomics analysis, see `omicverse-micro-metabol-paired`.

## Examples
- "Run the vsearch one-shot pipeline on the mothur MiSeq SOP FASTQs with `filter_max_ee=1.0` and the RDP v18 SINTAX DB; return an AnnData."
- "Re-run the same FASTQs with `backend='dada2'` and `filter_max_ee=2.0`; confirm diversity patterns match the vsearch result."
- "Compute Shannon + observed-OTUs alpha at the shallowest depth, then Bray-Curtis beta + 3-D PCoA, and plot the first two PCs colored by `obs['group']`."
- "Run Wilcoxon DA between Early and Late at phylum level with `min_prevalence=0.1` and report the top 10 hits."

## References
- Tutorial notebooks:
  - [`t_16s_amplicon.ipynb`](https://omicverse.readthedocs.io/en/latest/Tutorials-microbiome/t_16s_amplicon/) — vsearch end-to-end + alpha/beta/PCoA/DA on the mothur SOP.
  - [`t_16s_dada2.ipynb`](https://omicverse.readthedocs.io/en/latest/Tutorials-microbiome/t_16s_dada2/) — DADA2 backend on the same FASTQs.
- Live API verified — see [`references/source-grounding.md`](references/source-grounding.md).
