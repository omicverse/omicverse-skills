# Compatibility Notes

## Small-Dataset Failure

Current OmicVerse CytoTRACE2 has a real edge-case failure on chunks smaller than 100 cells.

Observed behavior:

- `process_subset(...)` skips KNN smoothing for chunks under 100 cells.
- the outer wrapper still expects a `CytoTRACE2_Score` column that is not produced on that branch
- the run then fails with `KeyError: 'CytoTRACE2_Score'`

Practical guidance:

- keep the working chunk size at 100 cells or more for smoke tests and production runs
- if a user truly has fewer than 100 cells, warn that the current wrapper may fail without an upstream fix

## Species Assumption

The wrapper uses gene-symbol casing heuristics to warn about species mismatches.

Practical guidance:

- mostly uppercase gene symbols generally indicate the human branch
- mixed-case mouse-style symbols generally indicate the mouse branch
- correct the species setting when the warning disagrees with the requested branch

## Parallelization Noise

Bounded CPU runs can emit repeated Numba or TBB fork warnings during CytoTRACE2 processing, especially when internal code paths spawn work under a threaded wrapper.

Practical guidance:

- prefer bounded CPU smoke runs with explicit `disable_parallelization=True`
- treat the warning as runtime noise unless it is followed by a hard failure

## Plotting Compatibility

The notebook uses `ov.utils.embedding(...)`, but the currently inspected interface is a compatibility wrapper with a generic signature.

Practical guidance:

- validate potency outputs first
- treat embedding overlay as optional presentation logic
