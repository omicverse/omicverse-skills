# Source Grounding

## Live Interfaces Checked

- `ov.single.DEG(adata, condition, ctrl_group, test_group, method='wilcoxon', use_raw=None)`
- `DEG.run(celltype_key, celltype_group=None, max_cells=100000, **kwargs)`
- `DEG.get_results()`

## Source-Verified Branches

### Constructor

- `method` supports `wilcoxon`, `t-test`, and `memento-de`.
- `use_raw=None` auto-detects `adata.raw`.
- If `use_raw` resolves to true and raw exists, the wrapper copies `adata.raw.X` and `adata.raw.var` into a fresh object before analysis.

### `run(...)`

- `celltype_group=None` expands to all unique values in the chosen cell-type column.
- The wrapper subsets to the selected cell types before any DEG backend runs.
- `max_cells` triggers downsampling only when the subset exceeds the threshold.
- For `wilcoxon` and `t-test`, the wrapper keeps only the control and test labels, then normalizes and log-transforms only when the matrix looks integer-valued.
- For `memento-de`, the wrapper creates a binary `stim` column, prefers the active matrix when it is integer-valued, otherwise falls back to a `counts` layer, and otherwise attempts count recovery from log-normalized values.
- Branch-specific kwargs are only forwarded to the memento call.

### `get_results()`

- The scanpy-backed branches return a standardized table with `log2FC`, `pvalue`, `padj`, `qvalue`, `size`, `sig`, `-log(pvalue)`, `-log(qvalue)`, `baseMean`, `pct_ctrl`, `pct_test`, and `pct_diff`.
- The memento branch returns the memento result table and then augments it with `baseMean`, `pct_ctrl`, `pct_test`, and `pct_diff`.

## Practical Branch Guidance

- Pick `wilcoxon` or `t-test` when the user wants the scanpy rank-gene path and is comfortable with the wrapper's normalization behavior.
- Pick `memento-de` only when the data contract is count-like enough to satisfy the branch's count assumptions.

## Reviewer-Side Evidence

- The notebook was read end to end and partitioned into DEG and DCT jobs before writing the skill.
- The installed OmicVerse source for `DEG` was inspected directly to confirm constructor defaults, `run(...)` branching, and `get_results()` postprocessing.
