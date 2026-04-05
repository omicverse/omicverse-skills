# Branch Selection

## Skill Boundary

The notebook starts with preprocessing, then fans out into multiple trajectory branches:

- preprocessing and PCA setup
- `TrajInfer(...).inference(method='diffusion_map')`
- `TrajInfer(...).inference(method='slingshot')`
- `TrajInfer(...).inference(method='palantir')`
- Palantir branch selection and gene trends
- `TrajInfer(...).inference(method='sctour')`

The preprocessing section is already a stable standalone job and belongs in the preprocessing skill.

The `diffusion_map`, `slingshot`, `palantir`, PAGA, and Palantir follow-ups stay together here because they all consume the same cluster-ready `AnnData`, share the same `TrajInfer` wrapper, and differ mostly by `method` plus a few follow-up calls.

`sctour` is split out because:

- it uses a different backend family
- it triggers a trainer-based workflow instead of a graph-only pseudotime branch
- it expects raw counts in `.X`
- it has a different compute and dependency profile

## `method='diffusion_map'`

Use this branch when:

- you want DPT-style pseudotime quickly
- you already trust the chosen `use_rep`
- you need a pseudotime key for downstream PAGA

This branch updates the AnnData in place and writes `dpt_pseudotime`.

## `method='slingshot'`

Use this branch when:

- you want explicit lineage curves fitted on an embedding
- the low-dimensional basis is biologically meaningful for curve fitting
- you can satisfy the extra `pcurvepy2` dependency

Notebook-specific note:

- `num_epochs=1` is only a short fitting path, not a universal default.
- `debug_axes` is optional and exists only for lineage visualization.

## `method='palantir'`

Use this branch when:

- you need entropy or fate probabilities, not just one pseudotime vector
- you know the starting state
- you may want lineage-specific branch masks or gene trends later

Prefer explicit `set_terminal_cells(...)` when you know the terminal states already.

Useful downstream branches:

- `palantir_cal_branch(...)` after `palantir_fate_probabilities` exists
- `palantir_cal_gene_trends(...)` after branch masks exist and an expression layer is ready

## PAGA After Trajectory Inference

PAGA is not a separate skill boundary here because in the notebook it is a topology-summary step over the same object.

Use it when:

- you want a cluster-level summary after `dpt_pseudotime`, `slingshot_pseudotime`, or `palantir_pseudotime`
- you need a coarse topology overlay for interpretation

Important parameter choices:

- `groups` should match the cluster key used in trajectory inference
- `use_time_prior` should point to the pseudotime key for the chosen branch
- `minimum_spanning_tree=True` keeps the tree-like summary path

## Palantir Gene Trends

Keep this inside the same skill instead of splitting it because:

- it is only useful after Palantir outputs exist
- it consumes Palantir branch masks and pseudotime
- it does not change the input contract, only the downstream result family

Do not present it as a required step for every Palantir run. It is optional and dependency-sensitive.
