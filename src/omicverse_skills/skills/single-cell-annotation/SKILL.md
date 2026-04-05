---
name: omicverse-single-cell-annotation
description: Annotate single-cell AnnData with OmicVerse using the CellTypist, gpt4celltype, or SCSA branches. Use when turning OmicVerse annotation notebooks into a reusable, triggerable skill, selecting a backend, or mapping a clustered AnnData object to cell-type labels.
---

# OmicVerse Single-Cell Annotation

## Goal

Annotate one clustered `AnnData` object with one OmicVerse backend at a time. Treat `celltypist`, `gpt4celltype`, and `scsa` as alternative branches that share the same input/output contract, not as sequential stages.

## Quick Workflow

1. If the input is raw counts, hand it off to the OmicVerse preprocessing path first; this skill starts from annotation-ready `AnnData` with a cluster column such as `leiden`.
2. Create `ov.single.Annotation(adata)` once.
3. Choose exactly one backend:
   - `celltypist` when you already have a CellTypist model path or want to resolve one with `query_reference(source="celltypist")`.
   - `gpt4celltype` when you want LLM-based cluster labels and can provide `AGI_API_KEY`.
   - `scsa` when you want marker-database scoring from clustered cells and have or can point to a local SCSA database.
4. Pass only the backend-specific kwargs for the selected branch.
5. Validate the expected `adata.obs` label column and any branch-specific side products before plotting or downstream analysis.

## Interface Summary

- `Annotation.annotate(method='celltypist', cluster_key='leiden', **kwargs)` dispatches across `celltypist`, `scsa`, `gpt4celltype`, `harmony`, `scVI`, and `scanorama`.
- `Annotation.query_reference(source='cellxgene', data_desc=None, llm_model='gpt-4o-mini', llm_api_key='sk*', llm_provider='openai', llm_base_url='https://api.openai.com/v1', llm_extra_params=None)` selects reference resources; source branches are `cellxgene` and `celltypist`.
- `Annotation.download_reference_pkl(reference_name, save_path, force_download=False)` and `Annotation.add_reference_pkl(reference)` support the CellTypist branch.
- `pySCSA(adata, foldchange=1.5, pvalue=0.05, output='temp/rna_anno.txt', model_path='', outfmt='txt', Gensymbol=True, species='Human', weight=100, tissue='All', target='cellmarker', celltype='normal', norefdb=False, cellrange=None, noprint=True, list_tissue=False, tissuename=None, speciename=None)` is the SCSA constructor.
- `pySCSA.cell_anno(clustertype='leiden', cluster='all', rank_rep=False)` and `pySCSA.cell_auto_anno(adata, clustertype='leiden', key='scsa_celltype')` expose the SCSA annotation flow.
- `gptcelltype(input, tissuename=None, speciename='human', provider='qwen', model='qwen-plus', topgenenumber=10, base_url=None)` is the LLM marker-to-label helper; provider branches are `openai`, `kimi`, and `qwen`.

Read `references/source-grounding.md` before documenting more detailed parameter behavior.

## Branch Selection

- Use `celltypist` when you want pretrained model annotation and already know the `.pkl` model, or when `query_reference(source="celltypist")` should choose one for you.
- Use `gpt4celltype` when you want prompt-driven annotation from cluster markers and can provide an API key.
- Use `scsa` when you want marker-database scoring from clustered cells and a local database path.
- Use `query_reference(source="celltypist")` only to choose a CellTypist model; do not treat it as a required stage.
- Treat `harmony`, `scVI`, and `scanorama` as source-present but out-of-scope extensions unless the user explicitly asks for transfer-learning annotation.

## Input Contract

- Start from `AnnData` that already has normalized expression and a cluster column such as `leiden`.
- For marker-based branches, keep `cluster_key` aligned with the label column you want to map.
- For CellTypist, you need a valid model path or a successful download step.
- For `gpt4celltype`, set `AGI_API_KEY` before the call; if you only want the generated prompt, call `gptcelltype(...)` directly.
- For SCSA, prefer a local database path instead of notebook-local absolute paths.
- If you still need raw preprocessing, run the notebook's QC, `preprocess`, `neighbors`, `leiden`, and `umap` handoff first.

## Minimal Execution Patterns

For CellTypist model selection and annotation:

```python
import omicverse as ov

obj = ov.single.Annotation(adata)
res = obj.query_reference(source="celltypist", data_desc="human PBMC scRNA-seq")
model_name = res.iloc[0]["model"]
obj.download_reference_pkl(model_name, save_path="models/celltypist.pkl")
obj.add_reference_pkl("models/celltypist.pkl")
obj.annotate(method="celltypist")
```

For LLM-based annotation:

```python
import os
import omicverse as ov

os.environ["AGI_API_KEY"] = "your-key"
obj = ov.single.Annotation(adata)
obj.annotate(
    method="gpt4celltype",
    cluster_key="leiden",
    tissuename="PBMC",
    speciename="human",
    model="gpt-5-mini",
    provider="openai",
    topgenenumber=5,
)
```

For prompt-only LLM smoke checks, call the helper directly:

```python
import os
from omicverse.single._gptcelltype import gptcelltype

os.environ["AGI_API_KEY"] = ""
prompt = gptcelltype({"0": ["MS4A1", "CD79A"]}, tissuename="PBMC", speciename="human")
```

For SCSA annotation:

```python
import omicverse as ov

obj = ov.single.Annotation(adata)
obj.add_reference_scsa_db("models/pySCSA.db")
obj.annotate(
    method="scsa",
    cluster_key="leiden",
    foldchange=1.5,
    pvalue=0.01,
    celltype="normal",
    target="cellmarker",
    tissue="All",
)
```

## Validation

After `celltypist`, check:

- `adata.obs["celltypist_prediction"]`
- `adata.obsm["celltypist_decision_matrix"]`
- `adata.obsm["celltypist_probability_matrix"]`

After `gpt4celltype`, check:

- `adata.obs["gpt4celltype_prediction"]` when the backend returned a real mapping
- the direct `gptcelltype(...)` helper returns a prompt string when you run the prompt-only path

After `scsa`, check:

- `adata.obs["scsa_prediction"]`

Before any branch, check:

- `cluster_key` exists in `adata.obs`
- the data object is already annotation-ready
- downstream plotting keys such as `X_umap` exist if you plan to render embeddings

## Constraints

- Do not use notebook-local absolute paths like `/scratch/...` in reusable docs or examples.
- Do not assume `celltypist` is installed in every environment.
- The prompt-only fallback lives in `gptcelltype(...)`, not in `Annotation.annotate(method="gpt4celltype")`.
- The live source exposes extra `method` branches `harmony`, `scVI`, and `scanorama`; they are not part of the notebook spine.
- For SCSA, keep the key you validate aligned with the branch you run: `pySCSA.cell_auto_anno(..., key='scsa_celltype')` uses its own default, while `Annotation.annotate(method='scsa', ...)` writes `scsa_prediction`.

## Resource Map

- Read `references/backend-selection.md` when choosing between the three annotation backends.
- Read `references/source-notebook-map.md` to trace notebook cells into this reusable skill.
- Read `references/source-grounding.md` for inspected signatures, source-level branch behavior, and compatibility notes.
- Read `references/compatibility.md` when notebook prose, parameter names, or current runtime behavior diverge.
