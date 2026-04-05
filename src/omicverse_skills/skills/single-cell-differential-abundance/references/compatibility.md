# Compatibility Notes

## Notebook Drift vs Current Wrapper Behavior

- The notebook notes that newer releases use a native Python implementation for the `milopy` branch while still keeping the legacy `milo` branch. That distinction matters because the two branches are still separate in current wrapper code.
- The notebook shows `run(num_samples=5000, num_warmup=500)` for scCODA. Current wrapper code forwards such kwargs only to the `sccoda` branch, not to `milopy` or `milo`.

## Runtime Considerations

- The Milo-family branches require a usable embedding in `obsm`; the notebook used Harmony, but any valid compatible embedding can satisfy that contract.
- `mix_threshold` is not a universal DCT option; it only changes neighborhood labels in the Milo-family result table.
- Keep acceptance and smoke execution shell-neutral. The generated skill only assumes a POSIX-like shell can invoke `${PYTHON} script.py`; it does not require `zsh` features.
