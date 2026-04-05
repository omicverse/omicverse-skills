# Branch Selection

## Constructor Choices

- Choose `feature_type='expression'` for the notebook-style expression workflow.
- Use a different `feature_type` only when the underlying modality really changes.
- Choose `ondisk=False` for bounded smoke runs and small datasets.
- Choose `ondisk=True` when large datasets force a disk-backed workflow.

## Post-Fit Choices

- Use `predicted(...)` for plain topic-derived labels.
- Use `get_results_rfc(...)` when the user wants classifier-derived labels from topic-confidence subsets.
- Choose `LDA_threshold` conservatively enough that each topic still has training examples.
