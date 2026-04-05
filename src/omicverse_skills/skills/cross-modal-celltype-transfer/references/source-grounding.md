# Source Grounding

## Inspected callables

- `omicverse.utils.weighted_knn_trainer(train_adata, train_adata_emb, n_neighbors=50)`
- `omicverse.utils.weighted_knn_transfer(query_adata, query_adata_emb, ref_adata_obs, label_keys, knn_model, threshold=1, pred_unknown=False, mode='package')`
- `omicverse.utils.mde(data, device=None, **kwargs)`
- `omicverse.utils.embedding(adata, basis, ..., show=None, save=None, return_fig=None, ...)`

## Source behavior worth preserving

- `weighted_knn_transfer` only implements `mode='package'`; any other value raises.
- `weighted_knn_transfer` matches `label_keys` by prefix against `ref_adata_obs.columns`, not by exact equality.
- `pred_unknown` controls whether `threshold` is used.
- `mde` imports `pymde` at runtime and raises if the package is missing.
- `embedding` can plot any existing basis in `.obsm`; it does not require MDE.

## Interface facts from inspection

- `weighted_knn_trainer` returns a `KNeighborsTransformer`.
- `weighted_knn_transfer` returns `(pred_labels, uncertainties)`.
- `mde` defaults to 2D output and chooses CPU or CUDA based on availability when possible.
