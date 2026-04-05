# Backend Selection

This skill keeps one boundary because the notebook's branches share one input contract and one output contract: a clustered `AnnData` object becomes cell-type labels in `adata.obs`. The branches are alternatives, not sequential jobs.

## Choose `celltypist`

Use when:

- you already have a pretrained CellTypist model
- you want to resolve a model name first with `query_reference(source="celltypist")`
- you want prediction matrices in `adata.obsm`

Notes:

- the branch needs the `celltypist` dependency
- the notebook's sample model is `Immune_All_Low.pkl`, but that is a worked example, not the identity of the skill

## Choose `gpt4celltype`

Use when:

- you want LLM-based cluster labeling from marker genes
- you can provide `AGI_API_KEY`
- you want the notebook's API-backed path or the prompt-only helper path

Notes:

- `Annotation.annotate(method="gpt4celltype")` assumes the helper can return a mapping
- if you only want the generated prompt, call `gptcelltype(...)` directly

## Choose `scsa`

Use when:

- you want marker-database scoring on clustered cells
- you already have a local SCSA database or can point to one
- you need `scsa_prediction` written into `adata.obs`

Notes:

- the constructor can accept `tissuename` and `speciename` aliases for compatibility
- the notebook's `/scratch/...` paths are not reusable and should be replaced with user-owned paths

## Skip The Other Method Branches

The live source also exposes `harmony`, `scVI`, and `scanorama` through `Annotation.annotate(...)`.

Treat them as a separate transfer-learning workflow unless the user explicitly asks for them. They were not part of the notebook's reusable spine.
