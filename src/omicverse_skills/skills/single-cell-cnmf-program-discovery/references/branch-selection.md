# Branch Selection

## Compute Mode

- Choose `use_gpu=False` for bounded smoke runs, CPU-only environments, or when GPU is unavailable.
- Choose `use_gpu=True` only for real accelerated runs in a compatible environment.

## Execution Mode

- Use `total_workers=1` for simple local runs.
- Use multiple workers only when you can run every worker index and later combine all outputs.

## Post-Consensus Choices

- Use `get_results(...)` for direct max-usage labels.
- Use `get_results_rfc(...)` when the user wants classifier-derived labels from thresholded usage programs.
