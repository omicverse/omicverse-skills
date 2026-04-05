# Compatibility Notes

## Compute And Filesystem Sensitivity

- cNMF can run in CPU or GPU mode; the notebook's broader conclusion about model quality does not change the fact that GPU availability is environment-specific.
- Persisted runs depend on a writable output location for intermediate and result files.
- Multi-worker runs are only complete if every worker index was executed.
- Prefer OmicVerse preprocessing wrappers such as `ov.pp.normalize_total`, `ov.pp.log1p`, and `ov.pp.pca` when preparing smoke or example inputs.

## Runtime Considerations

- Full cNMF runs are more expensive than simple graph clustering and should be smoke-tested on a reduced rank grid when validating locally.
- Keep acceptance and smoke execution shell-neutral.
