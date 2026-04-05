---
name: omicverse-single-cell-cytotrace2
description: Predict single-cell developmental potency with OmicVerse CytoTRACE2 on AnnData. Use when converting an OmicVerse CytoTRACE2 notebook into a reusable skill, when running the pretrained CytoTRACE2 potency workflow on mouse or human scRNA-seq data, or when choosing preprocessing mode, species, and parallelization settings for potency scoring and optional embedding overlays.
---

# OmicVerse Single-Cell CytoTRACE2

## Goal

Run the reusable CytoTRACE2 execution spine on single-cell `AnnData`: prepare the expression matrix, ensure pretrained weights are available, run `ov.single.cytotrace2(...)`, capture the score and potency columns written into `adata.obs`, and optionally overlay those outputs on an existing embedding. Keep the skill centered on the potency-prediction job rather than on the tutorial dataset.

## Quick Workflow

1. Inspect the input `AnnData`, especially whether `adata.X` contains the expression values you want CytoTRACE2 to score and whether the gene symbols match the declared species.
2. If the object still needs notebook-style normalization and HVG annotation, run `ov.pp.preprocess(...)` first and choose the `mode` branch explicitly.
3. Ensure the pretrained CytoTRACE2 weight directory is present before inference; download or stage it if needed.
4. Run `ov.single.cytotrace2(...)` with explicit `species`, batch settings, and parallelization settings.
5. Treat the returned table plus the five `adata.obs` columns as the shared handoff point for optional visualization.
6. Only if the user asks for plots, overlay `CytoTRACE2_Score`, `CytoTRACE2_Potency`, or `CytoTRACE2_Relative` on an existing embedding.
7. Validate the expected output columns and result table before treating the run as complete.

## Interface Summary

- `ov.pp.preprocess(adata, mode='shiftlog|pearson', target_sum=500000.0, n_HVGs=2000, organism='human', no_cc=False, batch_key=None, identify_robust=True)` adds a `counts` layer, robust-gene filtering, and HVG annotations.
- `mode` is a two-stage branch string split on `|`: the first token controls normalization (`shiftlog` or `pearson`), and the second controls HVG selection (`pearson` or `seurat`).
- `ov.single.cytotrace2(adata, use_model_dir, species='mouse', batch_size=10000, smooth_batch_size=1000, disable_parallelization=False, max_cores=None, max_pcs=200, seed=14, output_dir='cytotrace2_results')` runs pretrained CytoTRACE2 scoring on `adata.X`.
- `species` is a real branch with at least `mouse` and `human`; the implementation warns when the casing of gene symbols looks inconsistent with the declared species.
- `disable_parallelization` and `max_cores` control worker selection through `calculate_cores_to_use(...)`.
- `batch_size` controls chunking of the prediction stage; `smooth_batch_size` controls chunking of the KNN smoothing stage.
- Successful CytoTRACE2 runs write `CytoTRACE2_Score`, `CytoTRACE2_Potency`, `CytoTRACE2_Relative`, `preKNN_CytoTRACE2_Score`, and `preKNN_CytoTRACE2_Potency` into `adata.obs`, and also return a result table with those columns.
- The notebook uses `ov.utils.embedding(...)` for plotting; in current OmicVerse that is a compatibility wrapper, so treat embedding overlay as optional presentation rather than as part of the core potency job.

## Boundary

- Keep this as one skill because the notebook contains one tight job: optional preprocessing, CytoTRACE2 inference, then immediate inspection of the resulting potency columns.
- Do not split the embedding overlay into a separate skill unless the user already has CytoTRACE2 outputs and wants only a plotting helper.
- Do not absorb unrelated trajectory inference, annotation, or general preprocessing notebooks here.
- If a future notebook only changes `species`, batch sizes, or preprocessing `mode`, update this skill instead of creating a duplicate.

## Branch Selection

- Use `mode='shiftlog|pearson'` for the notebook path when you want log-normalization followed by Pearson-residual HVG selection.
- Use a first-token `pearson` normalization branch only when you explicitly want Pearson-residual normalization instead of shifted log normalization.
- Use a second-token `seurat` HVG branch only when you intentionally want Seurat-v3 HVG selection instead of Pearson-residual HVGs.
- Use `species='mouse'` for mixed-case mouse-style gene symbols and `species='human'` for mostly uppercase human-style symbols.
- Keep `disable_parallelization=False` for large runs when CPU parallelism is acceptable; set it to `True` for a simpler portable path or bounded smoke checks.
- Leave `max_cores=None` when automatic worker selection is acceptable; set it explicitly only when you need to cap concurrency.
- Keep `batch_size` at or below the dataset size; values above the cell count are reset internally.
- Avoid tiny datasets for production use. Current OmicVerse CytoTRACE2 fails on the branch where a chunk has fewer than 100 cells, because the skipped smoothing path does not produce the later `CytoTRACE2_Score` column expected by the wrapper.

## Input Contract

- Start from an `AnnData` whose `adata.X` is the matrix you want CytoTRACE2 to score.
- Ensure gene symbols are appropriate for the declared `species`.
- Keep an embedding such as `X_umap` only if you need the optional visualization stage; CytoTRACE2 itself does not require it.
- Ensure pretrained model weights are staged before the inference call.
- Expect intermediate files and a final tab-delimited result table to be written under the requested output directory.

## Minimal Execution Patterns

```python
import omicverse as ov

adata = ov.pp.preprocess(
    adata,
    mode="shiftlog|pearson",
    n_HVGs=2000,
)

results = ov.single.cytotrace2(
    adata,
    use_model_dir=model_dir,
    species="mouse",
    batch_size=10000,
    smooth_batch_size=1000,
    disable_parallelization=False,
    output_dir="cytotrace2_results",
)
```

```python
# Bounded smoke-style path
results = ov.single.cytotrace2(
    adata,
    use_model_dir=model_dir,
    species="mouse",
    batch_size=120,
    smooth_batch_size=60,
    disable_parallelization=True,
    max_pcs=50,
    output_dir="cytotrace2_smoke",
)
```

```python
# Optional embedding overlays after a successful run
ov.utils.embedding(
    adata,
    basis="X_umap",
    color=["CytoTRACE2_Score", "CytoTRACE2_Potency"],
    frameon="small",
)
```

## Validation

- Check that the model-weight directory exists before calling `ov.single.cytotrace2(...)`.
- Check that the returned table contains `CytoTRACE2_Score`, `CytoTRACE2_Potency`, `CytoTRACE2_Relative`, `preKNN_CytoTRACE2_Score`, and `preKNN_CytoTRACE2_Potency`.
- Check that those same five fields were written into `adata.obs`.
- Check that the output directory contains the final result table after the run.
- If you used a bounded smoke path, keep at least 100 cells so the current implementation reaches the working smoothing branch.
- If the declared `species` disagrees with the observed gene-symbol casing warnings, stop and correct the species assumption instead of trusting the scores.
- If you only validated the bounded CPU smoke path, say so explicitly and do not claim a full large-scale reproduction.

## Resource Map

- Read the branch-selection reference when choosing preprocessing `mode`, `species`, or concurrency settings.
- Read the source-grounding reference before extending the skill with more interface-specific claims.
- Read the notebook-map reference when deciding whether a future CytoTRACE2 tutorial belongs here or should become a separate downstream-only plotting skill.
- Read the compatibility reference when a run is small, parallelization behavior matters, or the wrapper fails on an edge case.
