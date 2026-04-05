# Integration Method Selection

## Why this is one skill

The notebook's stable job is cross-modal cell-type transfer from a reference AnnData to a query AnnData. The visual cells only inspect the same embeddings after transfer, so they should stay inside the same skill instead of becoming thin standalone skills.

## Core choice

- Use weighted KNN transfer as the core job.
- Treat `X_glue` as the notebook's example shared embedding, but keep the skill general to any shared basis already present in `.obsm`.
- Keep visualization optional.

## Branches

- `mode='package'` is the only runnable branch in the current source.
- `pred_unknown=True` enables threshold-based rejection to `"Unknown"`.
- `mde` is optional and depends on `pymde`; it is not part of the transfer contract.

## What to avoid

- Do not require the notebook's exact files or exact dataset.
- Do not require `X_umap` for the core path.
- Do not document `mode='paper'` as a usable branch unless the source changes.
