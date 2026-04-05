# Compatibility

## Extra Dependency

The `sctour` branch imports the external `sctour` package inside `TrajInfer.inference(method='sctour', ...)`.

Consequence:

- the branch cannot run unless that package is installed
- wrapper-level source grounding is still possible even when the runtime lacks the dependency

## Count Requirement

The wrapper initializes `sct.train.Trainer(..., loss_mode='nb', ...)`.

Consequence:

- raw UMI counts should be present in `.X`
- replacing `.X` with scaled or transformed values changes the intended loss assumptions

## Notebook-Specific Direction Flip

The notebook reverses pseudotime after training:

```text
adata.obs['sctour_pseudotime'] = 1 - adata.obs['sctour_pseudotime']
```

Treat that as a dataset-specific interpretation step, not as part of the universal contract.
