# Source Notebook Map

## Capability Partition

- Core executable job: fit an LDA topic model and write topic-derived labels.
- Optional downstream job: RFC-based hard labeling on top of fitted topic usage.
- Notebook-only pedagogy: long background on topic models and qualitative claims about one dataset.

## Why This Is A Separate Skill

- The notebook's graph clustering and cNMF sections are different model families with different runtime and dependency profiles.
- LDA topic modeling introduces MIRA-specific dependencies and an `ondisk` execution branch that do not belong in the graph clustering skill.

## What Was Demoted Out Of The Main Trigger Surface

- The dentate gyrus dataset.
- The notebook's threshold-picked topic count.
- The UMAP color panels for one run.
