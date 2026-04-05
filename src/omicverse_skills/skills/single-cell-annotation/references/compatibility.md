# Compatibility

Use this reference when notebook wording and current source behavior are not identical.

## Environment Constraints

- The `celltypist` dependency was not available in `omictest` during review.
- `gpt4celltype` can only produce real labels when the runtime has a usable `AGI_API_KEY`.
- Download-backed helpers such as `download_reference_pkl(...)` and `download_scsa_db(...)` may require network access.

## Notebook Versus Source

- The notebook sets `AGI_API_KEY` before the LLM branch; the direct helper path `gptcelltype(...)` is still the safer prompt-only smoke path.
- The notebook hardcodes `/scratch/...` and `temp/...` paths. Those are not reusable and should be replaced with user-owned paths.
- The notebook only demonstrates three `method` values, but live source also exposes `harmony`, `scVI`, and `scanorama`.

## Prompt-Only LLM Path

For prompt-only checks, set `AGI_API_KEY` to an empty string explicitly.

Do not rely on an unset environment variable to behave the same way, because the current helper checks for an empty string before deciding whether to short-circuit.

## SCSA Key Names

- `pySCSA.cell_auto_anno(..., key='scsa_celltype')` uses its own default output key.
- `Annotation.annotate(method='scsa', ...)` writes `scsa_prediction`.

Keep the validation key aligned with the call path you actually run.

## Conservative Rule

If notebook text and current source disagree:

1. trust the current source for executable behavior
2. mention the drift explicitly
3. keep the notebook's worked-example values out of the reusable trigger surface
