# Compatibility Notes

## Dependency And Runtime Sensitivity

- This skill depends on the MIRA stack being importable in the target environment.
- Current code applies a PyTorch compatibility shim for missing `torch._subclasses.schema_check_mode`.
- Full topic fitting can be materially slower than simple graph clustering and may be sensitive to PyTorch and accelerator configuration.
- Prefer OmicVerse preprocessing wrappers such as `ov.pp.normalize_total`, `ov.pp.log1p`, and `ov.pp.pca` when preparing smoke or example inputs.

## Filesystem Behavior

- `ondisk=True` writes local training and test datasets as part of execution.
- Keep acceptance and smoke execution shell-neutral.
