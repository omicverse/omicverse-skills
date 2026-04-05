---
name: omicverse-single-cell-scenic
description: Convert OmicVerse SCENIC notebooks into a reusable, triggerable skill for single-cell AnnData regulon analysis. Use when initializing SCENIC with cisTarget resources, choosing the RegDiffusion, GRNBoost2, or GENIE3 GRN branch, tuning regulon-construction thresholds, or running downstream RSS, binarization, and regulon-focused GRN exploration.
---

# OmicVerse Single-Cell SCENIC

## Goal

Run a reusable SCENIC analysis spine on single-cell `AnnData`: verify motif and ranking resources, initialize `ov.single.SCENIC`, infer a GRN, build regulons and AUCell scores, then optionally continue into RSS scoring, regulon binarization, and TF-centered GRN exploration. Keep this skill centered on one SCENIC result object rather than on the tutorial dataset.

## Quick Workflow

1. Inspect the input `AnnData`, especially the raw-count layer name, gene symbols, and grouping column for downstream interpretation.
2. Confirm that the ranking databases and motif annotation table exist before constructing `ov.single.SCENIC(...)`.
3. Pick the GRN branch up front with `method='regdiffusion'`, `method='grnboost2'`, or `method='genie3'`.
4. Run `SCENIC.cal_grn(...)` on a raw-count layer, then run `SCENIC.cal_regulons(...)` with explicit module-building kwargs when you need the notebook-style pruning behavior.
5. Treat the returned regulon `AnnData` plus `scenic_obj.auc_mtx`, `scenic_obj.regulons`, and `scenic_obj.modules` as the shared handoff point for all downstream stages.
6. Only if the user asks for interpretation, continue into RSS, binarization, embedding overlays, regulon marker ranking, or TF-target GRN plotting.
7. Validate the expected layers, tables, and object attributes before writing outputs.

## Interface Summary

- `ov.single.SCENIC(adata, db_glob, motif_path, n_jobs=8)` requires an `AnnData`, a ranking-database glob, and a motif annotation table.
- `SCENIC.cal_grn(method='regdiffusion'|'grnboost2'|'genie3', layer='counts', tf_names=None, **kwargs)` infers weighted TF-target edges and writes `edgelist` plus `adjacencies`.
- `SCENIC.cal_regulons(rho_mask_dropouts=True, seed=42, **kwargs)` builds modules, prunes them with cisTarget, computes AUCell scores, and returns a regulon `AnnData`.
- The forwarded module-building kwargs come from pySCENIC `modules_from_adjacencies(...)`, including `thresholds`, `top_n_targets`, `top_n_regulators`, `min_genes`, `absolute_thresholds`, `rho_dichotomize`, `keep_only_activating`, and `rho_threshold`.
- `regulon_specificity_scores(auc_mtx, cell_type_series)` converts the AUCell matrix into a cell-type-by-regulon RSS table.
- `binarize(auc_mtx, threshold_overides=None, seed=None, num_workers=1)` derives per-regulon thresholds and returns a binary activity matrix plus thresholds.
- `plot_rss(rss, cell_type, top_n=5, max_n=None, ax=None)` labels the top RSS regulons for one group.
- `sc.tl.rank_genes_groups(..., method='t-test'|'wilcoxon'|'logreg'|'t-test_overestim_var', use_raw=..., layer=...)` is the notebook's differential-ranking branch for regulon `AnnData`.
- `ov.single.cosg(..., mu=1, remove_lowly_expressed=False, expressed_pct=0.1, calculate_logfoldchanges=True, use_raw=True, layer=None)` is the faster marker-style alternative branch on the regulon matrix.
- `ov.single.build_correlation_network_umap_layout(...)`, `ov.single.add_tf_regulation(...)`, and `ov.single.plot_grn(...)` form the notebook's optional TF-target network visualization stage.

## Boundary

- Keep this as one SCENIC skill because the notebook's downstream stages all consume the same `SCENIC` outputs and are weak on their own without a completed regulon/AUCell result.
- Do not split just because a notebook changes `method`, threshold choices, dropout masking, or whether it continues into RSS or plotting.
- Split only when the user starts from an already computed regulon object or AUCell matrix and does not need SCENIC inference at all.
- Do not absorb unrelated preprocessing, trajectory inference, or cell-type annotation notebooks here.

## Branch Selection

- Use `method='regdiffusion'` when you want the notebook's neural GRN branch and can afford a training-oriented run.
- Use `method='grnboost2'` when you explicitly want the arboreto gradient-boosting branch, but verify the current wrapper on your runtime before treating it as production-ready.
- Use `method='genie3'` when you explicitly want that arboreto backend, and apply the same runtime verification caution as `grnboost2`.
- Use `layer='counts'` when the object follows the current live default; use `layer='raw_count'` only when the input object actually stores raw counts under that notebook-era name.
- Keep `rho_mask_dropouts=True` when reproducing the notebook's regulon path; turn it off only when you intentionally want the newer all-cells correlation behavior.
- Tune `thresholds`, `top_n_targets`, and `top_n_regulators` together because `modules_from_adjacencies(...)` combines all three module-construction strategies before pruning.
- Keep `keep_only_activating=True` for the portable path; set it to `False` only when you explicitly want both activating and repressing modules.
- Use `absolute_thresholds=False` unless you have a reason to express module thresholds in raw edge-weight units instead of percentiles.
- Use RSS when the user wants cell-type-specific regulons, binarization when they want on/off activity calls, and TF-network plotting only when they want a targeted exploratory view for selected TFs.
- Use `sc.tl.rank_genes_groups` when statistical testing matters; use `ov.single.cosg` when you want a fast specificity-style ranking branch on the regulon matrix.

## Input Contract

- Start from an `AnnData` with raw counts in a sparse layer that supports `.toarray()`.
- Ensure gene symbols are appropriate for the chosen motif/ranking resources.
- Provide a grouping column in `adata.obs` before RSS or marker-style downstream analysis.
- Expect `SCENIC.cal_grn(...)` to reject missing layers and to reject a layer that already looks log-normalized.
- Expect `SCENIC.cal_regulons(...)` to use `adata.to_df()` as the expression matrix for module pruning and AUCell scoring.
- Treat the motif table and ranking databases as external prerequisites, not as something this skill downloads.

## Minimal Execution Patterns

```python
import omicverse as ov

scenic = ov.single.SCENIC(
    adata,
    db_glob=db_glob,
    motif_path=motif_path,
    n_jobs=8,
)
edgelist = scenic.cal_grn(
    method="regdiffusion",
    layer="counts",
    n_steps=1,
    batch_size=8,
    device="cpu",
)
regulon_ad = scenic.cal_regulons(
    rho_mask_dropouts=True,
    thresholds=(0.75, 0.90),
    top_n_targets=(50,),
    top_n_regulators=(5, 10, 50),
)
```

```python
from omicverse.external.pyscenic.rss import regulon_specificity_scores
from omicverse.external.pyscenic.binarization import binarize

rss = regulon_specificity_scores(
    scenic.auc_mtx,
    scenic.adata.obs["cell_type"],
)
binary_mtx, auc_thresholds = binarize(
    scenic.auc_mtx,
    num_workers=1,
)
```

```python
import scanpy as sc
import omicverse as ov

sc.tl.rank_genes_groups(
    regulon_ad,
    groupby="cell_type",
    method="t-test",
    use_raw=False,
)
ov.single.cosg(
    regulon_ad,
    groupby="cell_type",
    key_added="cell_type_cosg",
    use_raw=False,
)
```

## Validation

- Check that the requested raw-count layer exists before `SCENIC.cal_grn(...)`.
- Check that `scenic.edgelist` and `scenic.adjacencies` are non-empty and contain `TF`, `target`, and `importance`.
- Check that `scenic.modules`, `scenic.regulons`, `scenic.auc_mtx`, and the returned regulon `AnnData` all exist after `SCENIC.cal_regulons(...)`.
- Check that the regulon `AnnData` obs index matches the source cells before copying embeddings or group labels onto it.
- Check that RSS output uses groups as rows and regulons as columns.
- Check that binarization returns both a binary matrix and threshold series with matching regulon names.
- Check that marker ranking wrote the expected key into `uns` before plotting.
- If you only validated the constructor plus a bounded `regdiffusion` smoke, say so explicitly and do not claim a full cisTarget reproduction.

## Resource Map

- Read the branch-selection reference when choosing the GRN backend, regulon module settings, or downstream interpretation branch.
- Read the source-grounding reference before extending the skill with more interface-specific claims.
- Read the notebook-map reference when deciding whether a new SCENIC tutorial belongs here or should become a separate downstream-only skill.
- Read the compatibility reference when a notebook uses `raw_count`, older dropout masking assumptions, or plotting helpers with extra dependencies.
