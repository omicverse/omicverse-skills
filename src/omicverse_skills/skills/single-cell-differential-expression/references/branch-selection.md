# Branch Selection

## Backend Choice

- Choose `wilcoxon` for the default fast nonparametric path.
- Choose `t-test` only when a t-test is requested explicitly or is statistically acceptable for the prepared matrix.
- Choose `memento-de` when the user wants memento specifically and the data contract can satisfy its count-oriented assumptions.

## Cell-Type Scope

- Choose one or a few explicit cell types when the question is lineage-specific.
- Choose all cell types only when the user really wants one pooled DEG pass rather than separate downstream comparisons.

## Matrix Source

- Let `use_raw=None` auto-detect raw only if you understand which matrix should drive the DEG run.
- Override `use_raw=False` when the active matrix is already the intended expression representation.
- For `memento-de`, prefer real counts or an explicit `counts` layer over relying on count recovery.

## Plotting

- Run violin plots only after the DEG table identifies genes worth checking.
- Treat plotting genes as request-specific output, not as part of the core reusable skill.
