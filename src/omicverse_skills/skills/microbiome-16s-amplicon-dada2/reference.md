# 16S amplicon — quick commands

## One-shot pipeline (vsearch backend)

```python
import omicverse as ov
import numpy as np, pandas as pd
ov.plot_set()

DB_FASTA = ov.alignment.fetch_rdp(db_dir='./db/rdp')

meta = pd.DataFrame({
    'group': [...],   # phenotype labels per sample
    'day':   [...],
}, index=[<sample_ids>])

adata = ov.alignment.amplicon_16s_pipeline(
    fastq_dir='./raw/MiSeq_SOP',
    workdir='./run_vsearch',
    db_fasta=DB_FASTA,
    threads=8, jobs=4,
    primer_fwd=None, primer_rev=None,    # if primers are NOT pre-trimmed, pass them
    filter_max_ee=1.0,
    unoise_minsize=2,
    sintax_cutoff=0.8,
    sample_metadata=meta,
)
print(adata)                       # samples × ASVs
```

## DADA2 backend — same wrapper

```python
adata_dada2 = ov.alignment.amplicon_16s_pipeline(
    samples=[
        ('F3D0', 'raw/F3D0_R1.fastq.gz', 'raw/F3D0_R2.fastq.gz'),
        # ...
    ],
    workdir='./run_dada2',
    db_fasta=DB_FASTA,
    backend='dada2',
    primer_fwd=None, primer_rev=None,
    filter_max_ee=2.0,
    sintax_cutoff=0.8,
    threads=4,
    sample_metadata=meta,
)
```

## Stepwise vsearch pipeline (when you want to inspect intermediates)

```python
# 1) merge paired reads
merge_res = ov.alignment.vsearch.merge_pairs(
    samples, output_dir='./step/merged',
    max_diffs=10, min_overlap=16,
    threads=8, jobs=4,
)

# 2) quality-filter (expected error)
filt_res = ov.alignment.vsearch.filter_quality(
    merge_res, output_dir='./step/filtered',
    max_ee=1.0, threads=8, jobs=4,
)

# 3) dereplicate → unique sequences
derep = ov.alignment.vsearch.dereplicate(
    filt_res, output_dir='./step/derep',
    min_uniq=2, threads=8,
)

# 4) UNOISE3 → ASVs
unoise = ov.alignment.vsearch.unoise3(
    derep['uniques'], output_dir='./step/asv',
    alpha=2.0, minsize=2, threads=8,
)

# 5) chimera removal
nochim = ov.alignment.vsearch.uchime3_denovo(
    unoise['asv'], output_dir='./step/asv',
)

# 6) SINTAX taxonomy
tax = ov.alignment.vsearch.sintax(
    nochim['asv'], db_fasta=DB_FASTA,
    output_dir='./step/taxonomy',
    cutoff=0.8, strand='both', threads=8,
)

# 7) per-sample OTU table
otutab = ov.alignment.vsearch.usearch_global(
    derep['combined'], nochim['asv'],
    output_dir='./step/otutab',
    identity=0.97, threads=8,
)

# 8) compose AnnData
adata = ov.alignment.build_amplicon_anndata(
    otutab_tsv=otutab['otutab'],
    asv_fasta=nochim['asv'],
    sintax_tsv=tax['tsv'],
    sample_metadata=meta,
    sample_order=[s[0] for s in samples],
)
```

## Alpha diversity

```python
min_depth = int(np.asarray(adata.X.sum(axis=1)).min())

ov.micro.Alpha(adata, rarefy_depth=min_depth).run(
    metrics=['shannon', 'observed_otus', 'simpson'],
)
print(adata.obs[['shannon', 'observed_otus', 'simpson']].describe())

# Single-metric convenience
shannon = ov.micro.Alpha(adata, rarefy_depth=min_depth).shannon()
observed = ov.micro.Alpha(adata, rarefy_depth=min_depth).observed()
```

## Beta diversity + PCoA / NMDS

```python
ov.micro.Beta(adata, rarefy_depth=min_depth).run(metric='braycurtis')
# distance matrix → adata.obsp['braycurtis']

ord_ = ov.micro.Ordinate(adata, dist_key='braycurtis')
ord_.pcoa(n=3)
# coords → adata.obsm['braycurtis_pcoa']
# var-fraction → adata.uns['micro']['braycurtis_pcoa_var']
pct = ord_.proportion_explained() * 100.0

# Or NMDS
ord_.nmds(n=2, random_state=0)
stress = adata.uns['micro']['braycurtis_nmds_stress']

# PCoA scatter
import matplotlib.pyplot as plt
coords = pd.DataFrame(adata.obsm['braycurtis_pcoa'],
                      index=adata.obs_names, columns=['PC1', 'PC2', 'PC3'])
fig, ax = plt.subplots(figsize=(5.5, 4.5))
for g, sub in adata.obs.groupby('group'):
    ax.scatter(coords.loc[sub.index, 'PC1'], coords.loc[sub.index, 'PC2'],
               label=g, s=70, alpha=0.85, edgecolor='k')
ax.set_xlabel(f'PC1 ({pct[0]:.1f}%)')
ax.set_ylabel(f'PC2 ({pct[1]:.1f}%)')
ax.legend(title='group'); plt.tight_layout(); plt.show()
```

## Differential abundance — Wilcoxon (default)

```python
da = ov.micro.DA(adata).wilcoxon(
    group_key='group', group_a='Early', group_b='Late',
    rank='phylum',                    # also: 'class', 'order', 'family', 'genus', 'species', or None for ASV
    min_prevalence=0.1,
)
print(da.head(10))                    # log2fc, pvalue, padj, mean_a, mean_b
```

## Differential abundance — pyDESeq2 (NB-GLM)

```python
da_deseq = ov.micro.DA(adata).deseq2(
    group_key='group', group_a='Early', group_b='Late',
    rank='genus',
    min_prevalence=0.1,
    alpha=0.05,
)
```

## Differential abundance — ANCOM-BC

```python
da_ancombc = ov.micro.DA(adata).ancombc(
    group_key='group',
    rank='genus',
    min_prevalence=0.1,
    pseudocount=1.0,
)
# (See the DA-comparison skill for trade-offs across the three methods.)
```

## Preprocessing helpers

```python
# Drop rare features
adata = ov.micro.filter_by_prevalence(adata, min_prevalence=0.1, min_count=1)

# Collapse ASVs to genus level
adata_genus = ov.micro.collapse_taxa(adata, rank='genus')

# CLR / ILR transforms (for Aitchison-style downstream)
adata = ov.micro.clr(adata, layer_out='clr')      # adata.layers['clr']
adata = ov.micro.ilr(adata, layer_out='ilr')

# Rarefy to a common depth (cache original counts at adata.layers['raw_counts'])
adata = ov.micro.rarefy(adata, depth=min_depth, seed=0,
                        save_original=True, copy=False)
```
