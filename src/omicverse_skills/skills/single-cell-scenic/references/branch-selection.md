# Branch Selection

## GRN Inference

- `SCENIC.cal_grn(method='regdiffusion', ...)` is the training-heavy neural branch.
- `SCENIC.cal_grn(method='grnboost2', ...)` is the lighter tree-based branch, but the current wrapper should be runtime-verified before you rely on it.
- `SCENIC.cal_grn(method='genie3', ...)` is the alternative arboreto tree branch and should be validated with the current wrapper before long runs.
- The live source default is `layer='counts'`. The tutorial uses `layer='raw_count'`, so choose the layer name from the actual object, not from notebook habit.
- `tf_names=None` leaves TF selection to the backend. Provide a TF list only when you intentionally want a restricted regulator set.

## Regulon Construction

- `SCENIC.cal_regulons(...)` forwards its extra kwargs into pySCENIC `modules_from_adjacencies(...)`.
- `thresholds` selects percentile or absolute-weight cutoffs for TF-target modules.
- `top_n_targets` adds top-target modules for each TF.
- `top_n_regulators` adds top-regulator-per-target modules.
- `min_genes` filters small modules after the TF itself is added.
- `absolute_thresholds=False` means `thresholds` are treated as percentiles, not raw edge weights.
- `rho_dichotomize=True` splits activating and repressing links using TF-target expression correlation.
- `keep_only_activating=True` keeps only activating modules on the portable path.
- `rho_threshold` controls when a TF-target correlation is considered activating or repressing.
- `rho_mask_dropouts=True` reproduces the tutorial path that ignores zero-expression cells when estimating TF-target correlations.

## Downstream Interpretation

- Use RSS when the question is “which regulons are specific to each cell group?”
- Use binarization when the question is “which cells are on versus off for each regulon?”
- Use embedding overlays when the question is “where is this regulon active on the manifold?”
- Use `sc.tl.rank_genes_groups` when statistical group-wise ranking matters and you need explicit `method` control.
- Use `ov.single.cosg` when you want a fast specificity-style ranking on the regulon matrix.
- Use TF-centered GRN plotting only after you already have selected TFs and a target dictionary from `scenic.regulons`.
