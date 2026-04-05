# Source Grounding

This skill is grounded in the current OmicVerse wrapper source that exposes the `sctour` branch.

## Primary Source File

- `omicverse/single/_traj.py`

## Inspected Interfaces

### `ov.single.TrajInfer`

Signature:

```text
(adata: anndata.AnnData, basis: str = 'X_umap', use_rep: str = 'X_pca', n_comps: int = 50, n_neighbors: int = 15, groupby: str = 'clusters')
```

### `TrajInfer.inference`

Signature:

```text
(self, method: str = 'palantir', **kwargs)
```

Observed source branches:

- `method='palantir'`
- `method='diffusion_map'`
- `method='slingshot'`
- `method='sctour'`

### `method='sctour'` Source Behavior

Observed wrapper logic in `omicverse/single/_traj.py`:

```text
import sctour as sct
tnode = sct.train.Trainer(self.adata, loss_mode='nb', **kwargs)
tnode.train()
self.adata.obs['sctour_pseudotime'] = tnode.get_time()
mix_zs, zs, pred_zs = tnode.get_latentsp(alpha_z=0.5, alpha_predz=0.5)
self.adata.obsm['X_TNODE'] = mix_zs
self.adata.obsm['X_VF'] = tnode.get_vector_field(self.adata.obs['sctour_pseudotime'].values, self.adata.obsm['X_TNODE'])
```

Grounded conclusions:

- the public branch selector is still `method='sctour'`
- the current wrapper hardcodes `loss_mode='nb'`
- the wrapper writes three notebook-relevant outputs: `sctour_pseudotime`, `X_TNODE`, and `X_VF`
- the notebook's `alpha_recon_lec` and `alpha_recon_lode` values are wrapper kwargs, not special syntax

## Runtime Finding

- the current environment could import OmicVerse but could not import the external `sctour` package itself

Consequence:

- this skill is source-grounded
- it was not fully empirically reproduced in the current local environment
