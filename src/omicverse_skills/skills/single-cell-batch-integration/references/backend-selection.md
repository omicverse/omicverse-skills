# Backend Selection Notes

## Default: Harmony

Choose `harmony` when:

- you already have PCA or are willing to let the helper derive it
- you want the notebook's main integration branch
- you want a representative CPU-friendly smoke path

Prefer `use_gpu=False` for deterministic lightweight validation unless GPU execution is explicitly required.

## Choose Combat

Choose `combat` when:

- you want a simple expression-space correction path
- you are fine with the helper rescoring and recomputing PCA internally
- you only need an integrated embedding written to `X_combat`

## Choose Scanorama

Choose `scanorama` when:

- the user explicitly asks for Scanorama
- the batch column is categorical and batch splitting is stable
- the environment has the extra Scanorama-side dependencies

## Choose scVI

Choose `scVI` when:

- raw counts are preserved in `adata.layers['counts']`
- the environment has `scvi-tools`
- you want a learned latent representation rather than a direct PCA correction
- you accept training cost that is materially higher than Harmony or Combat

Main notebook-aligned knobs:

- `n_layers`
- `n_latent`
- `gene_likelihood`

Live constructor also exposes:

- `dispersion`
- `latent_distribution`
- `dropout_rate`

## Choose CellANOVA

Choose `CellANOVA` when:

- you have a defensible `control_dict`
- you want the denoised reconstruction plus the derived PCA embedding
- you can provide enough control samples for the control pools to make sense

Main branch-specific knobs grounded from the underlying `calc_BE(...)` path:

- `control_dict`
- `var_cutoff`
- `k_max`
- `k_select`
- `verbose`

## Choose Concord

Choose `Concord` or `concord` when:

- Concord is installed
- you want its self-supervised integrated latent space
- you accept a heavier path than Harmony or Combat

Notebook-aligned kwargs such as `device='cuda'` are passed through the Concord constructor via `**kwargs`.

## Benchmarking Stage

Run benchmarking only when:

- you already have one or more integrated embedding keys
- you still have a non-integrated reference embedding such as `X_pca`
- both `batch_key` and `label_key` are present

For bounded smoke validation, reduce the metric bundle rather than trying to benchmark every possible metric on a large object.
