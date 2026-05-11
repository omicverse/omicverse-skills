---
name: bulk-fastq-quantification
title: Bulk RNA-seq FASTQ → count matrix with omicverse
description: End-to-end bulk RNA-seq quantification with omicverse's alignment module — SRA download, fastp QC, two interchangeable quantification paths (STAR + featureCount, OR alignment-free kb-python with technology='BULK'), and wiring into `ov.bulk.pyDEG` DESeq2. Single-cell kb-python (10XV2/10XV3) is out of scope — use the `single-cell-kb-alignment` skill instead.
---

## Overview

OmicVerse provides a complete bulk RNA-seq FASTQ-to-count-matrix pipeline via the `ov.alignment` module. This skill covers:

- **SRA data acquisition**: `prefetch` and `fqdump` (fasterq-dump wrapper), or `parallel_fastq_dump` (multi-threaded `.lite.1` → paired FASTQ)
- **Quality control**: `fastp` for adapter trimming and QC reports
- **Path A — alignment-based**: `STAR` aligner with auto-index building + `featureCount` (subread featureCounts wrapper) for gene quantification. Produces sorted BAMs and a gene × sample count DataFrame. Use when transcript structure / splicing matters or when downstream tools expect BAMs.
- **Path B — alignment-free (kb-python `technology='BULK'`)**: `ov.alignment.single.ref` builds a kallisto index + t2g + cDNA once per genome; `ov.alignment.count(technology='BULK', ...)` quantifies each sample directly to gene counts. Much faster than STAR + featureCounts on the same machine, with comparable gene-level accuracy for most DE use cases. This is the `t_mapping_kbpython` walkthrough.
- **DE handoff**: wiring per-sample outputs into `ov.bulk.pyDEG` (DESeq2 / edgeR / limma) → volcano plot.

All functions share a common CLI infrastructure (`_cli_utils.py`) that handles tool resolution, auto-installation via conda/mamba, parallel execution, and streaming output.

**Out of scope:** single-cell kb-python (`technology='10XV2'`, `'10XV3'`, velocity workflows). See `single-cell-kb-alignment` for that.

## Instructions

1. **Environment setup**
   - Bioinformatics tools are resolved automatically from PATH or the active conda environment.
   - If `auto_install=True` (default), missing tools are installed via mamba/conda on demand.
   - Supported tools for Path A: `prefetch`, `vdb-validate`, `fasterq-dump`, `fastp`, `STAR`, `samtools`, `featureCounts`, `pigz`, `gzip`.
   - For Path B (kb-python): `pip install kb-python`.

2. **SRA data download** (`ov.alignment.prefetch` + `ov.alignment.fqdump`, or `ov.alignment.parallel_fastq_dump`)
   - Use `prefetch` first for reliable downloads with integrity validation (`vdb-validate`).
   - Then convert to FASTQ with `fqdump`. It auto-detects single-end vs paired-end.
   - `fqdump` can also work directly from SRR accessions without prefetch.
   - For `.lite.1` direct downloads (e.g. NCBI mirrors), use `ov.datasets.download_data(url, dir=...)` + `ov.alignment.parallel_fastq_dump` instead.
   - All variants support retry with exponential backoff for network errors.
   ```python
   import omicverse as ov

   # Variant 1 — SRR accession via prefetch + fqdump
   pre = ov.alignment.prefetch(['SRR1234567', 'SRR1234568'],
                                output_dir='prefetch', jobs=4)
   fq  = ov.alignment.fqdump(['SRR1234567', 'SRR1234568'],
                              output_dir='fastq', sra_dir='prefetch',
                              gzip=True, threads=8, jobs=4)

   # Variant 2 — direct .lite.1 link + parallel_fastq_dump (paired-end)
   ov.datasets.download_data(
       'https://sra-downloadb.be-md.ncbi.nlm.nih.gov/.../SRR12544419.lite.1',
       dir='./data',
   )
   ov.alignment.parallel_fastq_dump(
       sra_id='./data/SRR12544419.lite.1',
       outdir='./data/SRR12544419',
       tmpdir='./tmp',
       threads=12, split_files=True, gzip=True,
   )
   ```

3. **FASTQ quality control** (`ov.alignment.fastp`)
   - Runs fastp for adapter trimming, quality filtering, and QC reporting.
   - Supports single-end and paired-end reads.
   - Produces per-sample JSON and HTML QC reports.
   - Sample format: tuple of `(sample_name, fq1_path, fq2_path_or_None)`.
   - Optional for Path B — kb-python can consume raw FASTQs directly, but fastp still improves DE quality on noisy SRA inputs.
   ```python
   samples = [
       ('S1', 'fastq/SRR1234567/SRR1234567_1.fastq.gz', 'fastq/SRR1234567/SRR1234567_2.fastq.gz'),
       ('S2', 'fastq/SRR1234568/SRR1234568_1.fastq.gz', 'fastq/SRR1234568/SRR1234568_2.fastq.gz'),
   ]
   clean = ov.alignment.fastp(samples, output_dir='fastp', threads=8, jobs=2)
   ```

---

### Path A — STAR + featureCounts (alignment-based)

4. **STAR alignment** (`ov.alignment.STAR`)
   - Aligns FASTQ reads using the STAR aligner.
   - **Auto-index building**: set `auto_index=True` (default) with `genome_fasta_files` and `gtf` to build index automatically if missing.
   - Produces coordinate-sorted BAM files.
   - Handles gzip-compressed FASTQs automatically (uses pigz/gzip/zcat).
   - Use `strict=False` (default) for graceful error handling per sample.
   ```python
   star_samples = [
       ('S1', 'fastp/S1/S1_clean_1.fastq.gz', 'fastp/S1/S1_clean_2.fastq.gz'),
       ('S2', 'fastp/S2/S2_clean_1.fastq.gz', 'fastp/S2/S2_clean_2.fastq.gz'),
   ]
   bams = ov.alignment.STAR(
       star_samples,
       genome_dir='star_index',
       output_dir='star_out',
       gtf='genes.gtf',
       genome_fasta_files=['genome.fa'],
       threads=8,
       memory='50G',
   )
   ```

5. **Gene quantification** (`ov.alignment.featureCount`)
   - Counts aligned reads per gene using featureCounts (subread).
   - Auto-detects paired-end from BAM headers (via pysam or samtools).
   - `auto_fix=True` (default) retries with corrected paired-end flag on error.
   - `gene_mapping=True` maps gene_id to gene_name from the GTF.
   - `merge_matrix=True` produces a combined count matrix across all samples.
   ```python
   bam_items = [
       ('S1', 'star_out/S1/Aligned.sortedByCoord.out.bam'),
       ('S2', 'star_out/S2/Aligned.sortedByCoord.out.bam'),
   ]
   counts = ov.alignment.featureCount(
       bam_items,
       gtf='genes.gtf',
       output_dir='counts',
       gene_mapping=True,
       merge_matrix=True,
       threads=8,
   )
   # counts is a pandas DataFrame (gene_id × samples) — feed directly to pyDEG
   ```

6. **Path A wiring**
   - **fastp → STAR**: fastp output is a list of dicts with keys `sample`, `clean1`, `clean2`, `json`, `html`.
     ```python
     star_samples = [
         (r['sample'], r['clean1'], r['clean2'] if r['clean2'] else None)
         for r in (clean if isinstance(clean, list) else [clean])
     ]
     ```
   - **STAR → featureCount**: STAR output is a list of dicts with keys `sample`, `bam` (or `error`).
     ```python
     bam_items = [
         (r['sample'], r['bam'])
         for r in (bams if isinstance(bams, list) else [bams])
         if 'bam' in r
     ]
     ```

---

### Path B — kb-python `technology='BULK'` (alignment-free)

7. **Build the kb reference** (`ov.alignment.single.ref`)
   - Builds the kallisto index, the transcript-to-gene map (`t2g.txt`), and the transcriptome FASTA (`cdna.fa`). Once per genome+annotation; reuse across the whole cohort.
   - Keep FASTA and GTF from the **same Ensembl release** to avoid transcript-ID mismatches.
   ```python
   ref_result = ov.alignment.single.ref(
       fasta_paths='genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz',
       gtf_paths='genomes/Homo_sapiens.GRCh38.108.gtf.gz',
       index_path='kb_ref/index.idx',
       t2g_path='kb_ref/t2g.txt',
       cdna_path='kb_ref/cdna.fa',
       temp_dir='tmp',
       overwrite=False,
   )
   ```

8. **Quantify each bulk sample** (`ov.alignment.count` with `technology='BULK'`)
   - `technology='BULK'` disables barcode/UMI handling — one observation per sample.
   - **Required bulk-specific args**: `filter_barcodes=False`, `parity='paired'|'single'`, `strand='unstranded'|'forward'|'reverse'`. The wrapper raises if `technology='BULK'` is combined with `filter_barcodes=True`.
   - Output per sample: `<output_path>/counts_unfiltered/adata.h5ad` (when `h5ad=True`) plus `cells_x_genes.genes.names.txt`.
   ```python
   for sra in ['SRR12544419', 'SRR12544421', 'SRR12544433', 'SRR12544435']:
       ov.alignment.count(
           fastq_paths=[
               f'./data/{sra}/{sra}.lite.1_1.fastq.gz',
               f'./data/{sra}/{sra}.lite.1_2.fastq.gz',
           ],
           index_path='kb_ref/index.idx',
           t2g_path='kb_ref/t2g.txt',
           technology='BULK',
           output_path=f'results/{sra}/',
           h5ad=True,
           filter_barcodes=False,
           parity='paired',
           strand='unstranded',
           threads=12,
       )
   ```

---

### Handoff to differential expression

9. **Wiring bulk outputs into `ov.bulk.pyDEG` (DESeq2)**

   Both paths end with a `gene × sample` integer count matrix that `pyDEG` consumes.

   **From Path A** — already a pandas DataFrame:
   ```python
   dds = ov.bulk.pyDEG(counts)         # counts from ov.alignment.featureCount(..., merge_matrix=True)
   ```

   **From Path B** — load per-sample h5ad, harmonise gene names, concatenate, transpose:
   ```python
   ad_dict = {}
   for sra in ['SRR12544419', 'SRR12544421', 'SRR12544433', 'SRR12544435']:
       ad = ov.read(f'./results/{sra}/counts_unfiltered/adata.h5ad')
       gene_name = ov.pd.read_csv(
           f'./results/{sra}/counts_unfiltered/cells_x_genes.genes.names.txt',
           header=None,
       )
       ad.var['gene_name'] = gene_name[0].tolist()
       ad.var['gene_id']   = ad.var.index
       ad.var.index        = ad.var['gene_name']
       ad.var_names_make_unique()
       ad.obs['sra']       = sra
       ad_dict[sra] = ad

   adata = ov.concat(ad_dict)
   adata.obs_names_make_unique()
   adata.obs['Group'] = ['no', 'no', 'yes', 'yes']   # phenotype labels

   counts = adata.to_df().T                          # gene × sample
   dds = ov.bulk.pyDEG(counts)
   ```

   **Run DE + filter + volcano** (same regardless of upstream path):
   ```python
   dds.drop_duplicates_index()
   result = dds.deg_analysis(
       treatment_groups=[...],   # sample IDs in the treatment arm (matrix columns)
       control_groups=[...],     # sample IDs in the control arm
       method='DEseq2',
   )
   # Optional: filter low-expression genes after DE
   result = result.loc[result['log2(BaseMean)'] > 1]

   dds.foldchange_set(fc_threshold=-1, pval_threshold=0.05, logp_max=10)
   dds.plot_volcano(title='DEG Analysis', figsize=(4, 4),
                    plot_genes_num=8, plot_genes_fontsize=12)
   ```

   - DESeq2 expects **raw integer counts** — do not log-normalize before `deg_analysis`.
   - `treatment_groups` / `control_groups` are sample IDs (matrix columns), not group labels.
   - For downstream GSEA / pathway / PPI, see the `gsea-enrichment`, `bulk-deg-analysis`, `bulk-stringdb-ppi` skills.

10. **Skipping completed steps**
    - All functions check for existing outputs and skip if `overwrite=False` (default).
    - Set `overwrite=True` to force re-execution.

11. **Troubleshooting**
    - If a tool is not found, check `auto_install=True` and that conda/mamba is accessible.
    - For STAR index errors, ensure `genome_fasta_files` points to uncompressed or gzip FASTA files.
    - For featureCounts paired-end detection errors, `auto_fix=True` handles most cases automatically.
    - GTF files can be gzip-compressed; they are auto-decompressed as needed.
    - For kb-python BULK: if quantification returns zero counts, double-check `parity` matches the actual library (single vs paired) and that `strand` matches the kit (most modern Illumina mRNA-seq is `unstranded` or `reverse`).

## Critical API Reference

### Sample format convention (Path A)

All alignment-based functions use a consistent sample tuple format:
- **FASTQ samples**: `(sample_name, fq1_path, fq2_path_or_None)`
- **BAM items**: `(sample_name, bam_path)` or `(sample_name, bam_path, is_paired_bool)`
- Single samples can be passed as a single tuple; multiple as a list of tuples.
- When a single tuple is passed, the return value is a single dict; for a list, a list of dicts.

### Auto-installation

```python
# All functions support these parameters:
auto_install=True   # Auto-install missing tools via conda/mamba
overwrite=False     # Skip if outputs already exist
threads=8           # Per-tool thread count
jobs=None           # Concurrent job count (auto-detected from CPU count)
```

## Examples

- **Bulk via STAR (alignment-based)**: `prefetch` → `fqdump` → `fastp` → `STAR` → `featureCount` → DataFrame → `ov.bulk.pyDEG` DESeq2 → volcano. (`t_mapping_STAR`)
- **Bulk via kb-python (alignment-free)**: `download_data` → `parallel_fastq_dump` → `ov.alignment.single.ref` → `ov.alignment.count(technology='BULK', filter_barcodes=False, parity='paired', strand='unstranded')` → per-sample h5ad → `ov.concat` → `ov.bulk.pyDEG` DESeq2 → volcano. (`t_mapping_kbpython`)
- **Local FASTQ + kb-python**: skip download, jump straight to `ov.alignment.single.ref` + `ov.alignment.count(technology='BULK', ...)`.
- **Local FASTQ + STAR**: skip download, start at `fastp` → `STAR` → `featureCount`.

## Related skills

- **`single-cell-kb-alignment`** — single-cell kb-python path (`technology='10XV2'` / `'10XV3'` / velocity workflows). Same `ov.alignment.single.ref` + `ov.alignment.count` API; different `technology` and barcode/UMI handling.
- **`bulk-deseq2-analysis`** — pure DESeq2 differential-expression analysis on an existing count matrix.
- **`bulk-deg-analysis`**, **`gsea-enrichment`**, **`bulk-stringdb-ppi`** — downstream DE/enrichment/network skills.

## References

- See [reference.md](reference.md) for copy-paste-ready code templates.
- Tutorial notebooks:
  - [`t_mapping_STAR.ipynb`](https://omicverse.readthedocs.io/en/latest/Tutorials-bulk/t_mapping_STAR/) — STAR + featureCounts walkthrough.
  - [`t_mapping_kbpython.ipynb`](https://omicverse.readthedocs.io/en/latest/Tutorials-bulk/t_mapping_kbpython/) — kb-python bulk walkthrough with DESeq2.
