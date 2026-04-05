# Source Notebook Map

Source notebook:

- `t_preprocess_cpu.ipynb`

Use this file to trace which notebook sections were preserved, condensed, or intentionally dropped.

## Notebook Sections To Skill Responsibilities

### 1. Imports And Runtime Setup

Notebook cells:

- `import scanpy as sc`
- `import omicverse as ov`
- `ov.plot_set(font_path='Arial')`
- `%load_ext autoreload`
- `%autoreload 2`
- `ov.settings.cpu_gpu_mixed_init()`

Preserved in the skill:

- the minimal OmicVerse import path
- the runtime toggle as an optional mode choice

Not preserved:

- auto-reload magics
- font-path setup as a hard requirement

### 2. Raw Data Load

Notebook cells:

- `sc.read_10x_mtx(...)`
- `adata.var_names_make_unique()`
- `adata.obs_names_make_unique()`

Preserved in the skill:

- the requirement that the caller provides an `AnnData` object
- the uniqueness step as a standard hygiene check

Not preserved:

- the notebook's 10x path and download instructions

### 3. QC

Notebook cells:

- `ov.pp.qc(...)`

Preserved in the skill:

- `ov.pp.qc(adata, **kwargs)` as the optional front-door for filtered input
- branch selection between `mode='seurat'` and `mode='mads'`
- doublet handling as a parameter, not as notebook narration

### 4. Preprocessing

Notebook cells:

- `ov.pp.preprocess(...)`
- `ov.pp.recover_counts(...)`
- raw-count inspection cells

Preserved in the skill:

- the core `preprocess -> scale -> pca -> neighbors` spine
- the `mode='shiftlog|pearson'` notebook path

Intentionally condensed or dropped:

- `recover_counts(...)` is treated as a diagnostic branch, not a required stage
- the notebook's demonstration of raw versus recovered CD3D values

### 5. Dimensionality Reduction

Notebook cells:

- `ov.pp.umap(...)`
- `ov.pp.mde(...)`
- `ov.pp.tsne(...)`
- `ov.pp.sude(...)`
- plotting calls for each embedding

Preserved in the skill:

- the fact that these are alternative embedding branches

Condensed:

- the skill does not force every plot to run in one fixed order
- use the embedding that matches the downstream question

### 6. Cell-Cycle Scoring

Notebook cells:

- `ov.pp.score_genes_cell_cycle(...)`
- `ov.pl.mde(adata_raw, color='phase')`

Preserved in the skill:

- a note that this is optional downstream analysis

Not elevated to a separate skill:

- it depends on the same preprocessed object, but it is not the main reusable job

### 7. Clustering And Labeling

Notebook cells:

- `ov.pp.leiden(...)`
- `ov.pl.embedding(..., color=['leiden', ...])`
- `ov.pl.ConvexHull(...)`
- `ov.utils.gen_mpl_labels(...)`

Preserved in the skill:

- `leiden` as the clustering stage before markers

Not preserved:

- the labeling and hull polish cells as mandatory steps

### 8. Marker Discovery

Notebook cells:

- `ov.single.find_markers(..., method='wilcoxon')`
- `ov.single.find_markers(..., method='cosg')`
- `ov.single.get_markers(...)`
- `ov.pl.markers_dotplot(...)`

Preserved in the skill:

- the marker branch with `wilcoxon` and `cosg`
- extraction with `get_markers`
- plotting with `markers_dotplot`

## Why This Stayed One Skill

The notebook contains multiple plot variants, but they all depend on the same preprocessing and clustering output. A separate plotting skill would repeat the same input contract and would not be independently useful until the object is already cluster-ready.

## What Was Intentionally Not Carried Over

- notebook-local absolute paths
- download commands
- debug-only cells
- figure-parity boilerplate
- example-specific output screenshots

## When To Reopen The Notebook

Reopen the notebook only when:

- the user wants exact figure parity
- the user wants the original 10x dataset rather than a synthetic smoke fixture
- a future OmicVerse release changes the branch names or storage keys
