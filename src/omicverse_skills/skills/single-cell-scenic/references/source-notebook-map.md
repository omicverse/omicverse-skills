# Source Notebook Map

## Capability Partition

### Core executable job

- Data readiness checks for the SCENIC resource files.
- `ov.single.SCENIC(...)` construction.
- GRN inference with `SCENIC.cal_grn(...)`.
- Regulon pruning and AUCell scoring with `SCENIC.cal_regulons(...)`.

### Optional downstream analysis jobs

- RSS calculation and RSS plotting.
- Regulon binarization and threshold visualization.
- Regulon-marker ranking with `sc.tl.rank_genes_groups` and `ov.single.cosg`.
- TF-centered GRN exploration using the regulon result plus additional graph helpers.

### Notebook-only material

- Dataset-specific plotting genes.
- Worked-example cell-type labels.
- One-off save and reload demonstration cells.

## Why This Stayed One Skill

- Every downstream block depends on the same completed `SCENIC` object, regulon list, and AUCell matrix.
- The notebook does not define a stable downstream-only input contract that is richer than “start from finished SCENIC outputs.”
- Splitting here would create a thin second skill whose main instructions would mostly be “load the existing SCENIC results and continue.”
- The boundary is still explicit: if a user already has regulon or AUCell outputs and does not need inference, that is the point where a future downstream-only skill would make sense.

## Section Mapping

- Notebook data-preparation and file-check cells map to the input-contract and compatibility notes.
- Notebook initialization and GRN cells map to the core execution spine and GRN branch rules.
- Notebook regulon and AUCell cells map to the regulon-construction branch rules and validation checks.
- Notebook RSS, binarization, and embedding cells map to optional downstream stage selection.
- Notebook TF-network cells map to the optional GRN-exploration stage, not to the minimum required workflow.
