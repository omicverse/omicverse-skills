# `ov.llm.SCLLMManager` — quick reference

## Constructor

```python
ov.llm.SCLLMManager(
    model_type: str = "scgpt",                # one of: scgpt | geneformer | scfoundation | uce | cellplm
    model_path: str | Path | None = None,     # checkpoint directory; required by every backend
    device: str | None = None,                # "auto" | "cpu" | "cuda" | "cuda:0"
    **kwargs,                                  # backend-specific (see below)
)
```

The constructor delegates to `ov.llm.ModelFactory.create_model(...)` and stores the wrapped backend on `self.model`. All instance methods route through `self.model`, so `model_type` decides what's supported.

### Backend-specific kwargs (most common)

- **scGPT** — `n_hvg: int = 1200`, `vocab_path: str | None = None`.
- **Geneformer** — `gene_id_col: str = "ensembl_id"`, `n_processes: int = 16`. *Input requires ENSEMBL gene IDs in `var`.*
- **UCE** — `species: str = "human"`. Set explicitly when running on non-human samples (zebrafish / mouse / macaque / pig / frog / mouse_lemur).
- **CellPLM** — `batch_size: int = 128` (largest of the five).
- **scFoundation** — `n_top_genes: int | None = None` (xTrimoGene HVG limit).

## Methods

### `get_embeddings(adata, **kwargs) -> np.ndarray`

Per-cell embedding; shape `(n_obs, n_embed_dim)`. Embedding dim depends on backend (scGPT 512, Geneformer 256, scFoundation 768, UCE 1280, CellPLM 512). Storage convention: assign to `adata.obsm["X_<backend>"]`.

### `annotate_cells(adata, **kwargs) -> dict`

Predicted cell-type labels. Returns `{"predictions": pd.Series, "probabilities": np.ndarray | None, "celltype_names": list[str]}`. Storage convention: write `predictions` into `adata.obs["cell_type_<backend>"]`.

- **scGPT / Geneformer** — require a fine-tuned annotation head (call `fine_tune(task="annotation")` first).
- **CellPLM** — zero-shot annotation supported via the pretrained head.
- **scFoundation** / **UCE** — annotation not exposed; raises `NotImplementedError`.

### `fine_tune(train_adata, valid_adata=None, task="annotation", **kwargs) -> dict`

Fine-tune the loaded backend on user data. `task` is one of `"annotation"` (cell-type classifier head) or `"integration"` (batch-aware encoder). Requires `train_adata.obs["celltype"]` for annotation, `train_adata.obs[batch_key]` for integration. Returns training-history dict (`{"loss": [...], "val_acc": [...], ...}`).

### `train_integration(train_adata, valid_adata=None, batch_key="batch", **kwargs) -> dict`

Convenience wrapper over `fine_tune(task="integration", batch_key=batch_key)`.

### `integrate(adata, batch_key="batch", **kwargs) -> np.ndarray`

Inference-only batch integration. Returns the integrated embedding `(n_obs, n_embed_dim)` and writes it to `adata.obsm["X_<backend>_integrated"]`. Use this when you've already trained / fine-tuned and want to apply the encoder to a new batch.

### `predict_celltypes(query_adata, **kwargs) -> dict`

Alias for `annotate_cells` on backends that expose a celltype head.

## Common errors

| Symptom | Cause | Fix |
|---|---|---|
| `Model 'scgpt' does not support task 'annotate'` (legacy) | You're calling the deprecated `ov.fm.run` validator. | Use `SCLLMManager(model_type='scgpt').annotate_cells(adata)` after `fine_tune(task='annotation')`. |
| `KeyError: 'ENSG...'` during Geneformer tokenization | Input genes are gene symbols, not ENSEMBL. | `ov.utils.geneid_to_ensembl(adata.var, ...)` before calling `get_embeddings`. |
| OOM in `get_embeddings` | Default batch_size too large for the GPU. | Pass `batch_size=8` (or lower) to `get_embeddings`. |
| `species not supported` (UCE) | Default `species='human'` mismatches data. | Pass `species='mouse'` (or whichever) at constructor time. |
| Empty `obsm` after `integrate` | Forgot to call `train_integration` first. | Train (or load a pre-trained integration head) before inference. |

## Patterns by use case

**"Embed cells with a foundation model and rebuild the kNN graph"**:
```python
manager = ov.llm.SCLLMManager(model_type="scgpt", model_path=PATH_SCGPT)
adata.obsm["X_scgpt"] = manager.get_embeddings(adata)
sc.pp.neighbors(adata, use_rep="X_scgpt")
sc.tl.umap(adata)
```

**"Zero-shot annotate PBMC with CellPLM"**:
```python
manager = ov.llm.SCLLMManager(model_type="cellplm", model_path=PATH_CELLPLM)
result = manager.annotate_cells(adata)
adata.obs["cell_type_cellplm"] = result["predictions"]
```

**"Compare two FM embeddings as obsm side-by-side"**:
```python
for tag, mt, path in [("scgpt", "scgpt", PATH_SCGPT), ("cellplm", "cellplm", PATH_CELLPLM)]:
    m = ov.llm.SCLLMManager(model_type=mt, model_path=path)
    adata.obsm[f"X_{tag}"] = m.get_embeddings(adata)
```
