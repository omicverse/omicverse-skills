# Compatibility

Use this reference when notebook wording and current source behavior diverge.

## Environment Constraints

- `harmony` is the most reliable representative smoke path and is light enough for a small synthetic test.
- `scVI` depends on `scvi-tools` and can be materially heavier than harmony.
- `scanorama` depends on extra integration packages and is source-present even though the notebook's main path is harmony.

## Notebook Versus Source

- The notebook uses a particular cellxgene reference dataset for convenience, but the reusable workflow should work with any query/reference pair that shares gene names.
- The notebook rewrites `adata_ref.var_names` from `feature_name`; that step is not universal, so only do it when your reference file actually needs it.
- The notebook demonstrates `shiftlog|pearson`, but the underlying preprocessing helper exposes additional mode combinations.

## Count Recovery Rule

If the query or reference already looks log-normalized, recover counts before concatenation.

The current constructor warns about this situation, so treating raw counts as the clean input contract is safer than assuming a preprocessed matrix will work unchanged.

## Key Naming Rule

Keep the validation key aligned with the branch you ran:

- `harmony_prediction` / `harmony_uncertainty`
- `scVI_prediction` / `scVI_uncertainty`
- `scanorama_prediction` / `scanorama_uncertainty`

The notebook also casts uncertainty scores to float before plotting. That is a plotting convenience, not a label-transfer requirement.

## Conservative Rule

If notebook text and current source disagree:

1. trust the current source for executable behavior
2. mention the drift explicitly
3. keep notebook-specific dataset names and file paths out of the reusable trigger surface
