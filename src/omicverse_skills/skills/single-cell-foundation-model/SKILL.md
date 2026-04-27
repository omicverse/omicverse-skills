---
name: omicverse-single-cell-foundation-model
title: Single-cell foundation model (SCLLMManager)
description: "Cell embedding, cell-type annotation, batch integration, and (where supported) perturbation prediction with single-cell foundation models — scGPT, Geneformer, scFoundation, UCE, CellPLM. Driven by the unified `ov.llm.SCLLMManager` interface; one object handles model loading, inference, and (optional) fine-tuning."
---

# Single-cell foundation model — `ov.llm.SCLLMManager`

## When to use this skill

Pick this skill when the user wants a transformer-style **per-cell representation** that is structurally richer than PCA: e.g. cell embedding, zero-shot or fine-tuned cell-type annotation, batch integration on the foundation-model latent, or (model permitting) perturbation prediction. The five fully-supported backends are **scGPT**, **Geneformer**, **scFoundation**, **UCE**, **CellPLM**; pick by data type, hardware, and gene-ID convention.

This skill **does not** cover marker-rule annotators (CellTypist / SCSA / gpt4celltype) or reference-mapping annotators (popV / scmap / SingleR). For those, prefer the `single-cell-annotation` skill — both can coexist; the foundation embedding is often used as the input space for downstream annotation/integration.

## Backend selection cheatsheet

| Model         | Tasks supported                          | Species                | Gene IDs | Min VRAM | CPU?  | Strengths                          |
|---------------|------------------------------------------|------------------------|----------|----------|-------|------------------------------------|
| **scGPT**     | embed, integrate, fine-tune→annotate     | human, mouse           | symbol   | 8 GB     | yes   | General RNA, longest-running, multi-modal extensions |
| **Geneformer**| embed, integrate, fine-tune→annotate     | human                  | ENSEMBL  | 4 GB     | **yes** | Ensembl-id pipelines, low-VRAM    |
| **scFoundation** | embed, integrate                       | human                  | symbol   | 16 GB    | no    | xTrimoGene architecture; perturbation work upstream |
| **UCE**       | embed, integrate                          | 7 species (cross-species) | symbol | 16 GB | no    | Zebrafish / macaque / pig / frog / lemur transfer |
| **CellPLM**   | embed, integrate, annotate (zero-shot)   | human                  | symbol   | 8 GB     | yes   | Fastest inference; cell-centric pretraining |

**Common pitfalls** the agent should avoid:

- *Don't* pick scGPT for an `annotate` task without first calling `fine_tune(task='annotation')` — scGPT's pretrained weights only emit embeddings; annotation is the fine-tuned head.
- *Don't* feed Geneformer gene **symbols** — its tokenizer expects ENSEMBL IDs. Convert via `ov.utils.geneid_to_ensembl(...)` first.
- *Don't* attempt `predict_perturbation` on scGPT/Geneformer/UCE/CellPLM — that path is scFoundation- or domain-specific. If perturbation is required, prefer the `single-cell-perturbation` skill.

## Canonical pipeline

```python
import omicverse as ov
import scanpy as sc

adata = sc.read_h5ad("query.h5ad")

# 1) Pick a backend + checkpoint path. SCLLMManager loads the model.
manager = ov.llm.SCLLMManager(
    model_type="scgpt",                 # or geneformer / scfoundation / uce / cellplm
    model_path="/path/to/checkpoints/scgpt",
    device="auto",                      # "cpu" | "cuda" | "auto"
)

# 2a) Cell embedding — works for all 5 backends.
embeddings = manager.get_embeddings(adata)
adata.obsm["X_scgpt"] = embeddings        # convention: X_<backend_short>

# 2b) Cell-type annotation — zero-shot (CellPLM) or via prior fine-tuning.
result = manager.annotate_cells(adata)
adata.obs["cell_type_scgpt"] = result["predictions"]

# 2c) Batch integration — uses the foundation latent as a batch-aware embedding.
manager.integrate(adata, batch_key="batch")
# adata.obsm["X_scgpt_integrated"] is populated.

# 2d) Optional fine-tune for annotation (when zero-shot accuracy isn't enough).
# train_adata must have obs["celltype"] ground-truth.
manager.fine_tune(train_adata, valid_adata, task="annotation")
result = manager.annotate_cells(adata)    # now uses fine-tuned head

adata.write_h5ad("query.h5ad")
```

## Storage conventions (for reviewers / graders)

After running the pipeline, the AnnData should carry:

- `obsm["X_<backend>"]` — the per-cell embedding (≥16-D, typically 256-1280-D depending on backend).
- `uns["neighbors"]` — kNN graph computed on `X_<backend>` (call `sc.pp.neighbors(adata, use_rep="X_<backend>")`).
- `obs["cell_type_<backend>"]` (annotation runs only) — predicted labels with the backend identifier in the column name so multi-method comparisons stay unambiguous.

## Discovery hooks

- `print(ov.utils.registry_lookup("scgpt"))` — surfaces `SCLLMManager` + its task-specific methods.
- `print(ov.utils.registry_lookup("foundation model embedding"))` — same.
- `print(ov.utils.registry_lookup("cell type annotation"))` does **not** intentionally surface SCLLMManager; that query routes to the rule-based / reference-mapping annotators in `single-cell-annotation`. Use the FM-explicit query when you mean the transformer path.

## References

- Tutorials (canonical examples): `Tutorials-llm/t_scgpt.ipynb`, `t_geneformer.ipynb`, `t_scfoundation.ipynb`, `t_uce.ipynb`, `t_cellplm.ipynb`.
- API: `omicverse.llm.SCLLMManager` and `omicverse.llm.ModelFactory` (low-level factory).
- Compatibility note (2026): the legacy `ov.fm.run` API was removed; SCLLMManager is the only supported foundation-model entry.

## Quick reference

See `reference.md` for the per-method signatures, return shapes, and caveats per backend.
