# Source Notebook Map

Source notebook:

- `t_anno_noref.ipynb`

Use this file to trace which notebook sections were preserved, condensed, or intentionally dropped.

## Notebook Sections To Skill Responsibilities

### 1. Import and plotting setup

Notebook cells:

- `import scanpy as sc`
- `import omicverse as ov`
- `ov.plot_set(font_path='Arial')`
- auto-reload magics

Preserved in the skill:

- not as a required step
- only the minimal OmicVerse imports that matter for execution

### 2. Raw-data preprocessing handoff

Notebook cells:

- `ov.pp.qc(..., tresh={...})`
- `ov.pp.preprocess(..., mode='shiftlog|pearson', n_HVGs=2000, target_sum=1e4)`
- `adata.raw = adata`
- `ov.pp.scale(adata)`
- `ov.pp.pca(...)`
- `ov.pp.neighbors(...)`
- `ov.pp.leiden(adata)`
- `ov.pp.umap(adata)`

Preserved in the skill:

- `SKILL.md` Input Contract
- `SKILL.md` Quick Workflow
- `references/compatibility.md`

Important source-grounded `mode` notes:

- `preprocess(mode='shiftlog|pearson')` is the notebook's selected handoff
- current source also exposes `shiftlog|seurat`, `pearson|pearson`, and `pearson|seurat`
- `qc(adata, **kwargs)` dispatches by runtime mode, with `seurat` and `mads` as the key QC branches

### 3. CellTypist branch

Notebook cells:

- `obj = ov.single.Annotation(adata)`
- `obj.query_reference(source='celltypist', ...)`
- `obj.download_reference_pkl(...)`
- `obj.add_reference_pkl(...)`
- `obj.annotate(method='celltypist')`
- `ov.pl.embedding(..., color='celltypist_prediction')`

Preserved in the skill:

- `SKILL.md` Interface Summary
- `SKILL.md` Minimal Execution Patterns
- `references/backend-selection.md`

### 4. gpt4celltype branch

Notebook cells:

- `os.environ['AGI_API_KEY'] = 'sk-*'`
- `obj.annotate(method='gpt4celltype', ...)`
- `ov.pl.embedding(..., color='gpt4celltype_prediction')`

Preserved in the skill:

- `SKILL.md` Minimal Execution Patterns
- `SKILL.md` Validation
- `references/compatibility.md`

### 5. SCSA branch

Notebook cells:

- `obj.download_scsa_db('temp/pySCSA_2024_v1_plus.db')`
- `obj.add_reference_scsa_db(...)`
- `obj.annotate(method='scsa', cluster_key='leiden', foldchange=1.5, pvalue=0.01, celltype='normal', target='cellmarker', tissue='All')`
- `ov.pl.embedding(..., color='scsa_prediction')`

Preserved in the skill:

- `SKILL.md` Branch Selection
- `SKILL.md` Validation
- `references/backend-selection.md`

## What Was Intentionally Not Carried Over

- notebook-local absolute paths
- debug-only cells such as `!pwd`
- repeated figures whose only purpose is notebook parity
- the notebook title as the skill identity

## When To Reopen The Notebook

Reopen the notebook only when:

- the user wants exact figure parity
- the user wants the original sample dataset rather than a reusable branch choice
- a future OmicVerse release changes the annotation branch contract
