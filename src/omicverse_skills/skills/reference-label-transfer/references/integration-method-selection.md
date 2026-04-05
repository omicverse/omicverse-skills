# Integration Method Selection

This skill covers one stable job: transfer labels from a reference AnnData to a query AnnData after building a shared integrated space.

## Choose `harmony`

Use when:

- you want the notebook's fast default transfer path
- you want the smallest dependency surface
- you want a representative method for smoke validation

Why this is the default:

- it is the first branch in both `AnnotationRef.train(...)` and `AnnotationRef.predict(...)`
- it produces the notebook's `harmony_prediction` and `harmony_uncertainty` outputs

## Choose `scVI`

Use when:

- you explicitly want the scVI latent space
- `scvi-tools` is installed
- you accept a heavier training path than harmony

Notes:

- the notebook demonstrates this branch explicitly
- the transfer keys become `scVI_prediction` and `scVI_uncertainty`

## Choose `scanorama`

Use when:

- you want Scanorama integration rather than Harmony or scVI
- the extra Scanorama dependencies are available
- you are okay with a branch that the notebook does not execute as its main path

Notes:

- the live source exposes the branch in both `train(...)` and `predict(...)`
- the transfer keys become `scanorama_prediction` and `scanorama_uncertainty`

## Preprocessing Mode

The notebook uses `mode='shiftlog|pearson'`.

The live preprocessing source also supports other mode combinations, so do not overfit the skill to the single notebook path:

- `shiftlog|seurat`
- `pearson|pearson`
- `pearson|seurat`

Use the notebook's mode only when reproducing that tutorial path exactly. Otherwise choose the mode that matches the data and downstream assumptions.
