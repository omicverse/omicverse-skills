---
name: omicverse-single-cell-kb-alignment
description: Build a triggerable kb reference and quantify single-cell FASTQs with OmicVerse alignment.single. Use when converting OmicVerse alignment notebooks into a reusable skill, when creating kallisto|bustools references with ref, or when running count on single-cell RNA-seq libraries.
---

# OmicVerse Single-Cell kb Alignment

## Goal

Build a kb reference from FASTA/GTF inputs, then quantify single-cell FASTQs into count matrices with OmicVerse alignment wrappers. This skill stops at the count matrix. The notebook's later QC, preprocessing, scaling, and PCA cells are generic follow-on work and should be handled separately.

## Quick Workflow

1. Choose the reference workflow and sequencing technology.
2. Build the reference with `ov.alignment.single.ref(...)`.
3. Run `ov.alignment.single.count(...)` on matching FASTQs.
4. Validate the returned paths and, if requested, load the emitted `adata.h5ad`.
5. Hand the count matrix to your downstream single-cell preprocessing workflow if needed.

## Interface Summary

- `ov.alignment.single.ref(index_path, t2g_path, fasta_paths=None, gtf_paths=None, cdna_path=None, workflow='standard', d=None, k=None, threads=8, overwrite=False, temp_dir='tmp', make_unique=False, include=None, exclude=None, dlist=None, dlist_overhang=1, aa=False, max_ec_size=None, nucleus=False, f2=None, c1=None, c2=None, flank=None, feature=None, no_mismatches=False, distinguish=False, **kwargs)` builds the reference bundle.
- `ov.alignment.single.count(index_path, t2g_path, technology, fastq_paths, output_path='.', whitelist_path=None, replacement_path=None, threads=8, memory='2G', workflow='standard', overwrite=False, temp_dir='tmp', tcc=False, mm=False, filter_barcodes=False, filter_threshold=None, loom=False, loom_names=None, h5ad=False, cellranger=False, gene_names=False, report=False, strand=None, parity=None, fragment_l=None, fragment_s=None, bootstraps=None, em=False, aa=False, genomebam=False, inleaved=False, batch_barcodes=False, exact_barcodes=False, numreads=None, store_num=False, long_read=False, threshold=0.8, platform='ONT', c1=None, c2=None, nucleus=False, **kwargs)` quantifies reads and returns discovered output paths.
- `ov.alignment.ref` and `ov.alignment.count` are aliases; `ov.alignment.single` is the convenient namespace used by the notebook.

## Boundary

- Keep this skill focused on kb reference building and count generation.
- Do not merge generic QC/preprocess/PCA into the same skill unless you want a separate downstream preprocessing skill.

## Branch Selection

- Use `workflow='standard'` unless the reference bundle requires a special mode.
- Use `nucleus=True` or `workflow='nucleus'` for nucleus/velocity workflows.
- Other `ref` workflows in the source include `lamanno`, `nac`, `kite`, and `custom`.
- Use `technology='10XV3'` or another preset that matches the library chemistry.
- `technology='BULK'` disables barcode filtering; `filter_barcodes` and `filter_threshold` must not be set there.
- In `count`, `tcc=True` and `mm=True` are mutually exclusive; the wrapper prefers `tcc`.
- Use `h5ad=True` when you want a quick AnnData artifact for downstream processing.

## Input Contract

- `ref` needs either `fasta_paths` plus `gtf_paths`, or a prebuilt reference bundle via `d`.
- `count` needs the reference `index_path` and `t2g_path` produced by `ref`.
- `fastq_paths` may be one file or a list of files; for paired-end 10x data, pass both R1 and R2 files.
- Keep `output_path` writable; the wrapper creates it if needed.

## Minimal Execution Patterns

Reference build:

```python
ref_info = ov.alignment.single.ref(
    index_path="ref/index.idx",
    t2g_path="ref/t2g.txt",
    fasta_paths="ref/transcripts.fa",
    gtf_paths="ref/genes.gtf",
)
```

Count matrix generation:

```python
count_info = ov.alignment.single.count(
    index_path="ref/index.idx",
    t2g_path="ref/t2g.txt",
    technology="10XV3",
    fastq_paths=["reads_R1.fastq.gz", "reads_R2.fastq.gz"],
    output_path="counts",
    h5ad=True,
    filter_barcodes=True,
)
```

## Validation

- `ref_info["index_path"]` and `ref_info["t2g_path"]` should exist after `ref`.
- `count_info` should include any output files that exist on disk, such as `h5ad_file`, `matrix_file`, `barcodes_file`, and `genes_file`.
- If `h5ad=True`, load the emitted `adata.h5ad` and confirm the expected shape.
- Confirm that `technology='BULK'` with `filter_barcodes=True` raises a `ValueError`.

## Resource Map

- Read `references/workflow-selection.md` for workflow and technology branch choices.
- Read `references/source-notebook-map.md` to see how the notebook maps to the reusable alignment boundary.
- Read `references/source-grounding.md` for inspected signatures and source behavior.
- Read `references/compatibility.md` for dependency and runtime caveats.
