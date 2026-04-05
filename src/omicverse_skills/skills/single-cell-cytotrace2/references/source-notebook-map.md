# Source Notebook Map

## Notebook To Skill Mapping

### Section: introduction and biological interpretation

- Kept out of the main skill body except for the minimal statement that CytoTRACE2 predicts developmental potency.
- The detailed potency-category explanation is notebook pedagogy, not required execution logic.

### Section: load dentate gyrus example

- Demoted to worked-example validation only.
- The reusable skill starts from a compatible `AnnData`, not from one specific tutorial dataset.

### Section: `ov.pp.preprocess(...)`

- Kept in the skill because the notebook executes it directly before CytoTRACE2 and because the `mode` branch materially changes behavior.

### Section: model download

- Kept in the skill as a prerequisite rule: pretrained weights must exist before inference.
- The exact hosting location is not the reusable skill identity, so only the need for staged weights remains in the main workflow.

### Section: `ov.single.cytotrace2(...)`

- Kept as the core job and documented from live source plus bounded execution evidence.

### Section: embedding overlays

- Kept as an optional post-inference branch.
- Not split into a second skill because it has no standalone input contract beyond the CytoTRACE2 output columns.

## Boundary Decision

This notebook does not justify multiple skills.

Reasons:

- preprocessing and CytoTRACE2 inference share one tight input object
- the plotting cells only consume the generated potency columns
- splitting would create a thin plotting-only wrapper with little independent value
