# Compatibility Notes

## Human Gene-Symbol Assumption

CellPhoneDB is built around human ligand-receptor resources.

Practical guidance:

- prefer human gene symbols
- verify the gene-symbol convention before treating CellPhoneDB output as biologically meaningful

## Small Iteration Trap

A reviewer-run reduced smoke path with `iterations=10` failed in the underlying CellPhoneDB code with `ZeroDivisionError`.

Practical guidance:

- do not use extremely small permutation counts for smoke tests
- keep the smoke path bounded with `threads=1` and a reduced cell subset, but still use a valid iteration setting such as `100`

## Input-Matrix Assumption

The tutorial explicitly checks that values look log-normalized and not scaled.

Practical guidance:

- avoid scaled expression matrices
- if the matrix looks inconsistent with the tutorial assumption, verify what CellPhoneDB expects for the specific analysis path before running

## Optional Plotting Dependencies

Some downstream views rely on optional plotting stacks beyond the core CellPhoneDB run.

Observed examples:

- Marsilea-backed bubble and heatmap views
- chord-diagram support for chord plots
- NetworkX-based centrality computations

Practical guidance:

- validate the core CellPhoneDB run first
- treat advanced plotting branches as optional extensions when runtime dependencies are available
