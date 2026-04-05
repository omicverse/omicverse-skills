---
name: omicverse-single-cell-lda-topic-clustering
description: Run OmicVerse single-cell LDA topic clustering with the MIRA backend as a reusable, triggerable skill. Use when fitting topic models on count-like AnnData, choosing between in-memory and on-disk execution, or converting topic assignments into hard cluster labels and RFC-based labels.
---

# OmicVerse Single-Cell LDA Topic Clustering

## Goal

Turn the notebook's MIRA-backed topic modeling section into a reusable job: fit an LDA topic model on count-like single-cell data, predict topic usage, and optionally derive hard cluster labels or RFC-based labels.

## Quick Workflow

1. Inspect whether the data still has a count layer and a highly-variable-feature flag.
2. Choose `feature_type`, `learning_rate`, and whether to run `ondisk=False` or `ondisk=True`.
3. Construct `LDA_topic(...)`, then decide how many topics to fit.
4. Run `predicted(num_topics=...)` to create `LDA_cluster`.
5. If the user wants harder classification on an embedding, run `get_results_rfc(...)` with explicit `use_rep` and `LDA_threshold`.

## Interface Summary

- `ov.utils.LDA_topic(adata, feature_type='expression', highly_variable_key='highly_variable_features', layers='counts', batch_key=None, learning_rate=0.001, ondisk=False)` constructs the wrapper.
- The documented `feature_type` branch is modality-sensitive; current docs describe at least `expression` and `accessibility`.
- `ondisk=False` keeps training data in memory. `ondisk=True` writes local training and test datasets before tuning and fitting.
- `predicted(num_topics=6)` fits the model with the chosen topic count, predicts topic usage, and writes `LDA_cluster`.
- `get_results_rfc(adata, use_rep='X_pca', LDA_threshold=0.5, num_topics=6)` trains tree-based classifiers on high-confidence topic assignments and writes `LDA_cluster_rfc` and `LDA_cluster_clf`.

## Stage Selection

- Use the plain `predicted(...)` path when the user wants one topic-to-cluster projection without extra classifiers.
- Use the RFC path only when the user explicitly wants classifier-derived hard labels from the topic representation.
- Use `ondisk=True` only when dataset size makes in-memory execution impractical.
- Treat the notebook's embedding panels as optional reporting, not the reusable contract.

## Input Contract

- Start from `AnnData`.
- Keep a count-like layer available under the selected `layers` key.
- Ensure the selected highly-variable flag exists in `adata.var`.
- Ensure the embedding named by `use_rep` exists before the RFC path.

## Minimal Execution Patterns

```python
import omicverse as ov

lda = ov.utils.LDA_topic(
    adata,
    feature_type="expression",
    highly_variable_key="highly_variable_features",
    layers="counts",
    batch_key=None,
    learning_rate=1e-3,
    ondisk=False,
)
lda.predicted(num_topics=13)
```

```python
lda.get_results_rfc(
    adata,
    use_rep="scaled|original|X_pca",
    LDA_threshold=0.4,
    num_topics=13,
)
```

## Constraints

- Do not use this skill on data that no longer has a usable count-like layer for the selected `layers` key.
- Do not treat the notebook's topic count as a universal default.
- Do not use the RFC path without a real embedding in `obsm`.
- `ondisk=True` creates local training/test datasets and is a runtime mode choice, not part of the reusable semantic contract.
- Keep smoke and acceptance commands shell-agnostic.

## Validation

- Check that the wrapper constructed successfully with the intended `feature_type`, `layers`, and `ondisk` mode.
- After `predicted(...)`, check that topic columns and `LDA_cluster` exist.
- After the RFC path, check that `LDA_cluster_rfc` and `LDA_cluster_clf` exist.
- If only a bounded smoke path was run, state whether model fitting was skipped or reduced.

## Resource Map

- Use the branch selection notes when choosing `feature_type`, `ondisk`, or the RFC path.
- Use the source grounding notes for current signatures, PyTorch compatibility behavior, and output labels.
- Use the notebook mapping notes to trace the notebook's topic-modeling section into this reusable skill.
- Use the compatibility notes for dependency, GPU, and filesystem-sensitive behavior.
