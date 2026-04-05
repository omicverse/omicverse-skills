# Compatibility Notes

- The wrappers depend on an external `kb` binary or an importable `kb` / `kb_python` module.
- The notebook's downstream QC/preprocess/PCA cells are not included here because they depend on the broader preprocessing stack and are not specific to kb alignment.
- `count` can emit `adata.h5ad`, `matrix.mtx`, `barcodes.tsv`, `genes.tsv`, and a `cellranger` directory when those outputs are requested and present.
- Keep `temp_dir` writable. The wrapper creates a unique temporary directory when the supplied path is unset, default, or already exists.
- For notebook-style paired-end data, pass the matching R1/R2 FASTQs together rather than trying to infer pairs inside the skill.
