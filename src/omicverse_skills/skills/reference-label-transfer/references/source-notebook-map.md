# Source Notebook Map

Source notebook:

- `t_anno_ref.ipynb`

Use this file to trace which notebook sections were preserved, condensed, or intentionally dropped.

## Notebook Sections To Skill Responsibilities

### 1. Reference-annotation setup

Notebook cells:

- `adata = ov.datasets.pbmc3k()`
- notebook prose explaining that reference annotation should use raw counts

Preserved in the skill:

- `SKILL.md` Input Contract
- `references/compatibility.md`

### 2. Reference selection and download

Notebook cells:

- `obj = ov.single.Annotation(adata)`
- `obj.query_reference(source='cellxgene', data_desc='PBMC for human', ...)`
- `ov.datasets.download_data_requests(...)`
- `adata_ref = ov.read(...)`
- `ov.pp.recover_counts(...)`
- `adata_ref.var.index = adata_ref.var['feature_name'].tolist()`

Preserved in the skill:

- `SKILL.md` Quick Workflow
- `SKILL.md` Input Contract
- `references/compatibility.md`

Important notebook-specific detail:

- the notebook downloads a particular cellxgene dataset for speed, but that dataset is a worked example, not the skill identity

### 3. Query/reference integration

Notebook cells:

- `objref = ov.single.AnnotationRef(adata_query=adata, adata_ref=adata_ref, celltype_key='cell_type')`
- `objref.preprocess(mode='shiftlog|pearson', n_HVGs=3000, batch_key='integrate_batch')`

Preserved in the skill:

- `SKILL.md` Interface Summary
- `SKILL.md` Branch Selection
- `references/integration-method-selection.md`

### 4. Harmony transfer branch

Notebook cells:

- `objref.train(method='harmony', n_pcs=50)`
- `ad_pre = objref.predict(method='harmony', n_neighbors=15)`
- `ov.pp.mde(ad_pre, use_rep='X_pca_harmony_anno')`
- `ov.pl.embedding(ad_pre, basis='X_mde', color='harmony_prediction')`
- `ad_pre.obs['harmony_uncertainty'] = ad_pre.obs['harmony_uncertainty'].astype(float)`

Preserved in the skill:

- `SKILL.md` Minimal Execution Patterns
- `SKILL.md` Validation

### 5. scVI transfer branch

Notebook cells:

- `objref.train(method='scVI', n_layers=2, n_latent=30, gene_likelihood='nb')`
- `ad_pre1 = objref.predict(method='scVI', n_neighbors=15)`
- `ov.pl.embedding(ad_pre1, basis='X_mde', color='scVI_prediction')`
- `ov.pl.embedding(ad_pre1, basis='X_mde', color='scVI_uncertainty', cmap='Reds')`

Preserved in the skill:

- `SKILL.md` Minimal Execution Patterns
- `SKILL.md` Validation
- `references/integration-method-selection.md`

## Source-Grounded Additions Beyond The Notebook

These details were added from live source and interface inspection, not just from notebook cells:

- `AnnotationRef` constructor signature and shared-gene check
- `AnnotationRef.preprocess(...)` signature and default batch key
- `AnnotationRef.train(...)` and `AnnotationRef.predict(...)` branch values
- the underlying `batch_correction(...)` helper exposes more methods than the notebook shows, including `combat`, `scanorama`, `CellANOVA`, and `Concord`

## What Was Intentionally Not Carried Over

- notebook-specific dataset download URLs
- notebook-only output file names
- notebook display prose
- tutorial-only cells that do not affect reusable execution

## When To Reopen The Notebook

Reopen the notebook only when:

- the user wants the exact cellxgene dataset from the tutorial
- the user wants notebook-parity figures instead of reusable transfer logic
- a future OmicVerse release changes the label-transfer branch names or output keys
