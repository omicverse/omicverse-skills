# Source Grounding

## Live Interfaces Checked

- `ov.utils.LDA_topic(adata, feature_type='expression', highly_variable_key='highly_variable_features', layers='counts', batch_key=None, learning_rate=0.001, ondisk=False)`
- `LDA_topic.predicted(self, num_topics=6)`
- `LDA_topic.get_results_rfc(self, adata, use_rep='X_pca', LDA_threshold=0.5, num_topics=6)`

## Source-Verified Branches

- `feature_type` is a modality branch; current docs describe `expression` and `accessibility`.
- The wrapper applies a PyTorch compatibility fix before importing MIRA.
- `ondisk=False` tunes and fits in memory.
- `ondisk=True` writes local train/test datasets and tunes against those on-disk datasets.
- `predicted(...)` fits the model, predicts topic usage, and writes `LDA_cluster` from the max topic per cell.
- `get_results_rfc(...)` trains decision-tree and random-forest classifiers from high-confidence topic subsets and writes `LDA_cluster_rfc` and `LDA_cluster_clf`.

## Practical Branch Guidance

- Use `ondisk=True` only when the data scale justifies the extra filesystem overhead.
- Use the RFC path only after topic prediction, not as a replacement for model fitting.

## Reviewer-Side Evidence

- The notebook was partitioned into graph clustering, LDA, cNMF, and evaluation before writing the skill.
- Live signatures and source behavior were checked for constructor, prediction, and RFC labeling methods.
