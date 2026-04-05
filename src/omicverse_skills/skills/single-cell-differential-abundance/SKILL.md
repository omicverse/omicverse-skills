---
name: omicverse-single-cell-differential-abundance
description: Run OmicVerse single-cell differential abundance or compositional analysis as a reusable, triggerable skill. Use when comparing cell-type abundance between conditions in AnnData, choosing between scCODA, milopy, and milo backends, or adapting a related notebook into a repeatable DCT workflow.
---

# OmicVerse Single-Cell Differential Abundance

## Goal

Turn notebook-style OmicVerse differential cell-type abundance analysis into a reusable execution spine. Keep this skill separate from DEG because DCT has a different constructor, different branch-specific prerequisites, and different result shapes.

## Quick Workflow

1. Inspect the input `AnnData`, the condition column, the cell-type column, and whether sample IDs and a usable embedding already exist.
2. Choose a backend before running anything: `sccoda` for the scCODA branch, `milopy` for the native Python Milo-like branch, or `milo` for the pertpy Milo branch.
3. For `milopy` and `milo`, make sure `sample_key` exists and `use_rep` points to an embedding already stored in `obsm`.
4. Construct `ov.single.DCT(...)` with explicit branch arguments, then run `run(...)`.
5. For `sccoda`, do any posterior-specific follow-up such as `set_fdr(...)` only after the model has finished.
6. Collect results with `get_results(...)` and verify the branch-specific output columns before plotting or exporting.

## Interface Summary

- `ov.single.DCT(adata, condition, ctrl_group, test_group, cell_type_key, method='sccoda', sample_key=None, use_rep=None)` initializes the wrapper.
- `method` branches are `sccoda`, `milopy`, and `milo`.
- `sample_key` is operationally required for the Milo-family branches because neighborhood counts are aggregated per sample.
- `use_rep` is required for `milopy`; current code also uses it for `milo`, even though only the `milopy` branch validates it eagerly.
- The constructor filters to the control and test labels immediately and reorders the condition categories.
- `DCT.run(**kwargs)` forwards kwargs only to the `sccoda` branch via `run_nuts(...)`.
- `milopy` and `milo` ignore extra kwargs in current wrapper code and instead run fixed neighborhood-testing logic.
- `DCT.get_results(mix_threshold=0.6)` returns a scCODA effect table for `sccoda`, or neighborhood-level Milo tables for `milopy` and `milo`.
- `mix_threshold` only affects the Milo-family branches by relabeling low-purity neighborhoods as `Mixed`.

## Stage Selection

- Use `method='sccoda'` when the user wants compositional Bayesian inference and is prepared for posterior sampling and FDR tuning.
- Use `method='milopy'` when the user wants the native Python Milo-like path and already has a valid embedding such as PCA or Harmony.
- Use `method='milo'` only when the pertpy Milo branch is specifically required.
- Keep diagnostic plots and beeswarm visualizations optional; they are downstream reporting on top of the DCT results, not the reusable contract itself.

## Input Contract

- Start from an `AnnData` object with condition labels and cell-type annotations in `obs`.
- Provide control and test labels that exist in the condition column.
- For `sccoda`, provide a sample identifier because the scCODA preparation step loads cell-level data into a sample-aware modality.
- For `milopy` and `milo`, provide both `sample_key` and `use_rep`.
- Ensure the embedding named by `use_rep` already exists in `adata.obsm`.

## Minimal Execution Patterns

```python
import omicverse as ov

dct = ov.single.DCT(
    adata,
    condition="condition",
    ctrl_group="Control",
    test_group="Salmonella",
    cell_type_key="cell_label",
    method="sccoda",
    sample_key="batch",
)
dct.run(num_samples=5000, num_warmup=500)
res = dct.get_results()
```

```python
dct = ov.single.DCT(
    adata,
    condition="condition",
    ctrl_group="Control",
    test_group="Salmonella",
    cell_type_key="cell_label",
    method="milopy",
    sample_key="batch",
    use_rep="X_harmony",
)
dct.run()
res = dct.get_results(mix_threshold=0.6)
```

## Constraints

- Do not assume one abundance backend can replace another without changing the input contract.
- Do not hide `method`, `sample_key`, or `use_rep`; those choices control materially different code paths.
- Do not pass `num_samples` or `num_warmup` as if they affect the Milo-family branches; current wrapper code only forwards kwargs to `sccoda`.
- Do not claim that neighborhood purity relabeling matters for `sccoda`; `mix_threshold` only changes the Milo-family result table.
- Keep smoke and acceptance commands shell-agnostic. Do not rely on `zsh`-specific syntax, shell startup files, or shell-only features when a direct `${PYTHON} script.py` command is enough.
- Do not put local paths or environment names into the reusable instructions.

## Validation

- Check that the control and test labels both remain after constructor filtering.
- For `sccoda`, check that `get_results()` returns an effect table and that any FDR adjustment is applied after model fitting.
- For `milopy` or `milo`, check that the result table includes neighborhood annotation columns and abundance statistics such as FDR-like outputs.
- Check that the embedding named by `use_rep` existed before construction.
- If only one representative branch was run as a smoke test, say which branch was validated and which branches were source-grounded only.

## Resource Map

- Use the branch selection notes when choosing between scCODA and the Milo-family branches.
- Use the source grounding notes when you need live signatures, current branching behavior, or branch-specific caveats.
- Use the notebook map when tracing which tutorial cells became which reusable instruction.
- Use the compatibility notes when notebook prose and current wrapper behavior differ.
