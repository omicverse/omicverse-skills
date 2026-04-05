# Source Grounding

## Live Interfaces Checked

- `ov.single.DCT(adata, condition, ctrl_group, test_group, cell_type_key, method='sccoda', sample_key=None, use_rep=None)`
- `DCT.run(**kwargs)`
- `DCT.get_results(mix_threshold=0.6)`

## Source-Verified Branches

### Constructor

- `method` supports `sccoda`, `milopy`, and `milo`.
- The constructor filters to the selected control and test labels immediately.
- `sccoda` builds a pertpy scCODA object and prepares a cell-level modality keyed by the condition.
- `milopy` requires `use_rep` explicitly and validates that the embedding key exists.
- `milo` also uses `use_rep` during neighbor construction, even though the constructor does not validate it as strictly as `milopy`.

### `run(...)`

- `sccoda` forwards kwargs to `run_nuts(...)`, then calls `credible_effects(...)`.
- `milopy` runs fixed neighborhood testing, graph building, and annotation without consuming extra kwargs.
- `milo` also runs fixed neighborhood testing, graph building, and annotation without consuming extra kwargs.

### `get_results(...)`

- `sccoda` returns the scCODA effect table directly.
- `milopy` and `milo` both relabel neighborhoods whose top annotation fraction falls below `mix_threshold` as `Mixed`.
- `mix_threshold` therefore matters only for the Milo-family branches.

## Practical Branch Guidance

- Choose `sccoda` when posterior sampling and effect inclusion are part of the request.
- Choose `milopy` when the user wants the OmicVerse native Python branch and already has an embedding.
- Choose `milo` only when the pertpy Milo backend is required explicitly.

## Reviewer-Side Evidence

- The notebook was read end to end and partitioned into DEG and DCT jobs before writing the skill.
- The installed OmicVerse source for `DCT` was inspected directly to confirm constructor defaults, backend branching, and the branch-specific meaning of `run(...)` kwargs and `mix_threshold`.
