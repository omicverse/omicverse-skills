---
name: omicverse-single-cell-differential-expression
description: Run OmicVerse single-cell differential expression analysis as a reusable, triggerable skill. Use when comparing conditions inside one or more cell types in AnnData, choosing between Wilcoxon, t-test, and memento backends, or adapting a related notebook into a repeatable DEG workflow.
---

# OmicVerse Single-Cell Differential Expression

## Goal

Turn notebook-style OmicVerse DEG analysis into a compact execution spine for condition-vs-condition comparisons inside selected cell types. Keep this skill focused on DEG only; differential abundance and compositional analysis are a separate job with a different input contract and backend surface.

## Quick Workflow

1. Inspect the input `AnnData`, the condition column, the cell-type column, and whether raw counts are available in `adata.raw` or `adata.layers["counts"]`.
2. Choose a DEG backend before running anything: `wilcoxon` or `t-test` for the scanpy rank-gene path, `memento-de` for the count-aware memento path.
3. Subset to the required cell types with `celltype_key` and `celltype_group`; leave `celltype_group=None` only when you really want all cell types pooled into one DEG run.
4. Run `DEG.run(...)` with an explicit `max_cells` policy and method-specific kwargs when using memento.
5. Collect results with `get_results()`, then verify the expected statistical columns and expression-percentage columns before saving or plotting.

## Interface Summary

- `ov.single.DEG(adata, condition, ctrl_group, test_group, method='wilcoxon', use_raw=None)` initializes the wrapper.
- `method` branches are `wilcoxon`, `t-test`, and `memento-de`.
- `use_raw=None` auto-detects `adata.raw`; if raw exists, the wrapper copies `adata.raw.X` into a fresh object and uses that matrix for DEG.
- `DEG.run(celltype_key, celltype_group=None, max_cells=100000, **kwargs)` subsets by cell type and dispatches to the selected backend.
- `wilcoxon` and `t-test` normalize and log-transform only when the active matrix still looks integer-valued.
- `memento-de` expects count-like input; if `X` is not integer-valued, the wrapper first looks for `adata.layers["counts"]` and otherwise tries count recovery from log-normalized values.
- memento-specific kwargs are passed through `run(...)`; the notebook used `capture_rate`, `num_cpus`, and `num_boot`.
- `DEG.get_results()` returns a DEG table and augments it with `baseMean`, `pct_ctrl`, `pct_test`, and `pct_diff`.

## Stage Selection

- Use `method='wilcoxon'` when you want the notebook's fast nonparametric path on one or more selected cell types.
- Use `method='t-test'` only when a t-test is acceptable for the data distribution and normalization state.
- Use `method='memento-de'` when the user explicitly wants memento's differential-mean path or needs the extra branch-specific kwargs.
- Keep `celltype_group` explicit for targeted comparisons such as one lineage or cluster; set it to `None` only when a pooled all-cell-type DEG run is actually intended.
- Treat violin plots as optional reporting after DEG, not as part of the reusable execution contract.

## Input Contract

- Start from an `AnnData` object with a condition column in `obs`.
- Provide control and test labels that exist in that condition column.
- Provide a cell-type column for `celltype_key`.
- Prefer preserving raw counts in `adata.raw` or `adata.layers["counts"]`, especially for `memento-de`.
- Expect the wrapper to coerce dense matrices to CSR internally.

## Minimal Execution Patterns

```python
import omicverse as ov

deg = ov.single.DEG(
    adata,
    condition="condition",
    ctrl_group="Control",
    test_group="Salmonella",
    method="wilcoxon",
)
deg.run(
    celltype_key="cell_label",
    celltype_group=["TA"],
    max_cells=100000,
)
res = deg.get_results()
```

```python
deg = ov.single.DEG(
    adata,
    condition="condition",
    ctrl_group="Control",
    test_group="Salmonella",
    method="memento-de",
    use_raw=False,
)
deg.run(
    celltype_key="cell_label",
    celltype_group=["TA"],
    capture_rate=0.07,
    num_cpus=12,
    num_boot=5000,
)
res = deg.get_results()
```

## Constraints

- Do not merge DEG and DCT into one skill at execution time; they are independently triggerable and use different backends.
- Do not assume the notebook's example cell type or example genes are universal defaults.
- Do not leave `method` implicit; branch choice changes both the statistical engine and the data expectations.
- Do not claim that `memento-de` ran on raw counts unless you actually confirmed raw counts, a `counts` layer, or successful count recovery.
- Keep smoke and acceptance commands shell-agnostic. Do not rely on `zsh`-specific syntax, shell startup files, or shell-only features when a direct `${PYTHON} script.py` command is enough.
- Do not put local paths or environment names into the reusable instructions.

## Validation

- Check that the selected control and test labels both remain after subsetting.
- Check that `DEG.run(...)` used the intended cell-type subset.
- For `wilcoxon` and `t-test`, check that the result table contains at least `log2FC`, `pvalue`, `padj`, `qvalue`, `baseMean`, `pct_ctrl`, `pct_test`, and `pct_diff`.
- For `memento-de`, check that the result table still gains `baseMean`, `pct_ctrl`, `pct_test`, and `pct_diff` after `get_results()`.
- If the dataset is large, say whether `max_cells` caused downsampling.
- If only a bounded smoke path was run, say so explicitly.

## Resource Map

- Use the branch selection notes when choosing the DEG backend or deciding whether all cell types should be pooled.
- Use the source grounding notes when you need the live signatures, defaults, or branch-specific behavior.
- Use the notebook map when tracing which tutorial cells became which reusable instruction.
- Use the compatibility notes when notebook prose and current runtime behavior diverge.
