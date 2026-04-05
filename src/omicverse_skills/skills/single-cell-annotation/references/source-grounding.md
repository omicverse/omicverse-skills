# Source Grounding

This skill is grounded in the live OmicVerse source and inspected runtime signatures.

## Inspection Notes

- Source inspected from the local OmicVerse checkout.
- Interface inspection ran in `omictest`.
- The annotation workflow was checked against live signatures and source branches before being written into the skill.

## Primary Source Files

- `omicverse/single/_annotation.py`
- `omicverse/single/_gptcelltype.py`
- `omicverse/single/_anno.py`
- `omicverse/pp/_preprocess.py`
- `omicverse/pp/_qc.py`

## Inspected Signatures

### `Annotation`

Inspected signature:

```text
(adata: anndata.AnnData)
```

This class stores one query `AnnData` object and exposes annotation helpers.

### `Annotation.annotate`

Inspected signature:

```text
(self, method='celltypist', cluster_key='leiden', **kwargs)
```

Source-grounded branch values:

- `celltypist`
- `scsa`
- `gpt4celltype`
- `harmony`
- `scVI`
- `scanorama`

Notebook-covered branches:

- `celltypist`
- `gpt4celltype`
- `scsa`

### `Annotation.query_reference`

Inspected signature:

```text
(self, source='cellxgene', data_desc: str = None, llm_model='gpt-4o-mini', llm_api_key='sk*', llm_provider='openai', llm_base_url='https://api.openai.com/v1', llm_extra_params=None)
```

Source-grounded branch values:

- `source`: `cellxgene`, `celltypist`
- `llm_provider`: `openai`, `custom_openai`, `doubao`, `anthropic`, `ollama`

### `Annotation.download_reference_pkl`

Inspected signature:

```text
(self, reference_name: str, save_path: str, force_download: bool = False) -> str
```

### `Annotation.add_reference_pkl`

Inspected signature:

```text
(self, reference: str)
```

### `Annotation.add_reference_scsa_db`

Inspected signature:

```text
(self, reference: str)
```

### `Annotation.download_scsa_db`

Inspected signature:

```text
(self, save_path: str = 'temp/pySCSA_2023_v2_plus.db')
```

### `gptcelltype`

Inspected signature:

```text
(input, tissuename=None, speciename='human', provider='qwen', model='qwen-plus', topgenenumber=10, base_url=None)
```

Source-grounded provider branches:

- `openai`
- `kimi`
- `qwen`

### `pySCSA`

Inspected signature:

```text
(adata: anndata.AnnData, foldchange: float = 1.5, pvalue: float = 0.05, output: str = 'temp/rna_anno.txt', model_path: str = '', outfmt: str = 'txt', Gensymbol: bool = True, species: str = 'Human', weight: int = 100, tissue: str = 'All', target: str = 'cellmarker', celltype: str = 'normal', norefdb: bool = False, cellrange: str = None, noprint: bool = True, list_tissue: bool = False, tissuename: str = None, speciename: str = None) -> None
```

Inspected methods:

```text
cell_anno(self, clustertype: str = 'leiden', cluster: str = 'all', rank_rep=False) -> pandas.DataFrame
cell_auto_anno(self, adata: anndata.AnnData, clustertype: str = 'leiden', key='scsa_celltype') -> None
```

## Source-Derived Behavior Notes

### CellTypist branch

- `Annotation.annotate(method='celltypist')` calls `celltypist.annotate(...)` with `majority_voting=True`
- it writes `celltypist_prediction`, `celltypist_decision_matrix`, and `celltypist_probability_matrix`
- `download_reference_pkl(...)` resolves a model name to a downloadable URL and stores the model path back on the object

### gpt4celltype branch

- `Annotation.annotate(method='gpt4celltype')` calls `get_celltype_marker(...)` and then `gptcelltype(...)`
- the notebook's `provider='openai'` and `model='gpt-5-mini'` values are valid on the helper side
- prompt-only mode is exposed by `gptcelltype(...)` when `AGI_API_KEY` is set to an empty string

### SCSA branch

- `Annotation.annotate(method='scsa')` uses `cluster_key='leiden'` by default and writes `scsa_prediction`
- `pySCSA` also supports the lower-level `cell_anno(...)` and `cell_auto_anno(...)` flow
- `pySCSA` supports compatibility aliases `tissuename` and `speciename`

### Preprocessing handoff

The notebook's annotation input is produced by the OmicVerse preprocessing path:

- `preprocess(adata, mode='shiftlog|pearson', target_sum=500000.0, n_HVGs=2000, organism='human', no_cc=False, batch_key=None, identify_robust=True)`

The live source also exposes `mode` combinations that the notebook does not display explicitly:

- `shiftlog|seurat`
- `pearson|pearson`
- `pearson|seurat`

The QC helper dispatches by runtime mode as well, with CPU, mixed, and GPU paths.

## Interface Inspection Evidence

Key branch terms were confirmed in the live source and by inspection:

- `method` for `Annotation.annotate`
- `source` and `llm_provider` for `Annotation.query_reference`
- `provider` for `gptcelltype`
- `mode` for `preprocess`

## Human Review Notes

- The notebook's CellTypist model download path is not reusable as written because it hardcodes a local filesystem location.
- The `celltypist` dependency was not available in the `omictest` runtime used for review, so the CellTypist branch could not be fully smoke-tested there.
