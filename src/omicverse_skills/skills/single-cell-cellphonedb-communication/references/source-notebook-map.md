# Source Notebook Map

## Notebook To Skill Mapping

### Section: CellPhoneDB background and motivation

- Demoted out of the main skill body.
- Useful for notebook reading, but not part of the reusable execution contract.

### Section: load annotated example data

- Demoted to worked-example validation.
- The reusable skill starts from a compatible annotated `AnnData`, not from one fixed example file.

### Section: run `ov.single.run_cellphonedb_v5(...)`

- Kept as the core job.
- Documented from live source and bounded execution rather than from notebook prose alone.

### Section: persist and reload results

- Demoted from the main skill.
- Persistence is optional workflow glue, not the stable analysis identity.

### Section: initialize `CellChatViz`

- Kept because every downstream notebook branch depends on the same visualization-ready interaction object and optional palette handling.

### Section: aggregated circle plots

- Kept as an optional high-level summary branch.

### Section: pathway-level aggregation and significance

- Kept as an optional downstream branch because it still consumes the same `CellChatViz` object and is a realistic user request on top of CellPhoneDB results.

### Section: ligand-receptor contribution, bubble, chord, and signaling-role analysis

- Kept inside the same skill as optional downstream branches.
- Not split because they share the same input contract and object state.

## Boundary Decision

One skill is better than multiple small skills here.

Reasons:

- the notebook revolves around one processed CellPhoneDB result object
- every downstream branch depends on the same `adata_cpdb` layout and interaction metadata
- splitting by visualization type would mostly duplicate the same trigger surface and setup steps
