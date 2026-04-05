# Source Notebook Map

## Capability Partition

- Core executable job: differential cell-type abundance or compositional analysis between two conditions.
- Optional downstream job: DCT-specific diagnostic plots and beeswarm-style visualizations after abundance testing.
- Notebook-only material: the worked example dataset and the long explanatory prose around one intestinal infection scenario.

## Why This Is A Separate Skill

- The same notebook also contains DEG analysis, but DEG is independently triggerable and uses a different wrapper with different statistical assumptions.
- DCT introduces a distinct `method` surface and branch-specific prerequisites such as `sample_key` and `use_rep`.

## Notebook To Skill Mapping

- The scCODA section became one explicit backend branch with posterior-sampling notes.
- The milopy section became a separate backend branch with embedding and sample requirements.
- The legacy milo section stayed as an optional backend branch because the wrapper still exposes it.
- Plotting cells were demoted to optional reporting because the reusable contract is the abundance test itself.

## Example-Specific Material That Was Removed From The Trigger Surface

- The Haber dataset loader and exact disease labels.
- The specific color palette for beeswarm plots.
- The notebook's hardcoded examples of one cell-type fraction threshold for demonstration.
