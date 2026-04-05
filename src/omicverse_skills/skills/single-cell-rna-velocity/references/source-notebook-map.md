# Source Notebook Map

## Capability Partition

- Core reusable job: single-cell RNA velocity analysis on a velocity-ready `AnnData`
- Comparison job: run the same velocity pipeline across `scvelo`, `dynamo`, and `latentvelo`
- Refinement job: apply `graphvelo` on top of an existing velocity layer and project the refined vectors
- Display-only job: the final UMAP embedding and stream plot

## Cell Map

- Cells 0-2: imports and dataset loading
- Cells 4-11: `latentvelo` branch followed by `graphvelo` refinement
- Cells 13-16: `scvelo` branch followed by `graphvelo` refinement
- Cells 18-24: `dynamo` branch followed by `graphvelo` refinement

## Boundary Decision

- The notebook is a downstream GraphVelo refinement tutorial on top of earlier velocity outputs, not a new analysis family.
- That makes it a better fit for updating the existing RNA velocity skill than for creating a duplicate skill.
- The skill boundary stays at "velocity-ready AnnData in, prior velocity layer refined and projected out".
- The final stream plot remains in scope because it is the user-visible completion of each branch.

## Handoff

- If you need raw read quantification, use `omicverse-single-cell-kb-alignment` first.
- If you already have a velocity-ready `AnnData`, continue with this skill.
