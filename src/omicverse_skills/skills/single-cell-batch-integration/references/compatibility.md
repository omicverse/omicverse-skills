# Compatibility Notes

- `harmony` and `combat` are the best representative smoke paths for a bounded validation run.
- `scanorama` depends on extra packages such as `intervaltree` and `fbpca`.
- `scVI` depends on `scvi-tools` and requires `adata.layers['counts']`.
- `CellANOVA` needs a valid `control_dict` and a preprocessed object that is already normalized, log-transformed, and scaled.
- `Concord` is available in this environment, but it is materially heavier than the Harmony and Combat smoke paths.
- `ov.pp.mde(...)` and `ov.utils.mde(...)` are optional visualization helpers, not proof that integration succeeded.
- `Benchmarker.prepare()` computes neighbors at multiple neighborhood sizes; keep smoke datasets small.
- For bounded acceptance, it is better to run a reduced benchmark metric bundle on a tiny synthetic dataset than to pretend the full notebook benchmark was validated.
- The reusable skill should stay shell-agnostic and should not assume any specific shell startup file or shell-only syntax.
