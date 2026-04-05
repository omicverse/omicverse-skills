# Source Notebook Map

## Capability Partition

- Core executable job: run one clustering backend on a prepared single-cell embedding or graph.
- Optional downstream jobs: embedding plots and ARI comparison across backend outputs.
- Notebook-only pedagogy: long clustering background and discussion of which method "looks best" on one dataset.

## Why This Is One Skill

- Leiden, Louvain, scICE, and GMM all answer the same user request: choose a clustering backend for a prepared `AnnData`.
- They share one tight input contract around embeddings and graph availability.
- The LDA and cNMF sections do not belong here because they are different model families with different runtime profiles and outputs.

## What Was Demoted Out Of The Main Trigger Surface

- The dentate gyrus dataset identity.
- The notebook's exact UMAP panels.
- The final ARI table comparing notebook-specific labels.
