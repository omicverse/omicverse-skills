# Source Notebook Map

## Notebook cells to reusable steps

- Import and plot settings: setup only.
- `wget` FASTA/GTF downloads: notebook-specific data acquisition, not part of the reusable skill.
- `ov.alignment.single.ref(...)`: reusable reference build.
- `ov.alignment.single.count(...)`: reusable quantification step.
- `sc.read_h5ad(...)` plus gene-name assignment: reusable validation of count output.
- `ov.pp.qc(...)`, `ov.pp.preprocess(...)`, `ov.pp.scale(...)`, `ov.pp.pca(...)`, and the final `ov.pl.embedding(...)`: generic downstream preprocessing that is intentionally outside this skill boundary.

## Why the split

The notebook contains one alignment/counting job and one generic preprocessing job. The alignment job is the reusable notebook-specific capability; the later preprocessing cells can be applied to any AnnData produced elsewhere, so they belong in a separate downstream workflow rather than in this skill.
