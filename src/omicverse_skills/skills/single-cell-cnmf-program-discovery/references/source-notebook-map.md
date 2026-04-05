# Source Notebook Map

## Capability Partition

- Core executable job: fit consensus NMF programs and load normalized usage results.
- Optional downstream job: RFC-based hard labeling from usage programs.
- Notebook-only pedagogy: qualitative comparison to the other clustering methods and dataset-specific figure panels.

## Why This Is A Separate Skill

- cNMF has a distinct compute profile, an explicit CPU/GPU branch, and a multi-stage factorization workflow.
- It is independently triggerable from graph clustering and LDA topic modeling.

## What Was Demoted Out Of The Main Trigger Surface

- The dentate gyrus example.
- The notebook's specific `selected_K` and `density_threshold`.
- The notebook's program-usage panels.
