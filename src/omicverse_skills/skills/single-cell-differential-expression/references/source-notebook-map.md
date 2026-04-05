# Source Notebook Map

## Capability Partition

- Core executable job: condition-vs-condition differential expression within selected cell types.
- Optional downstream job: violin plots for selected marker genes after DEG.
- Notebook-only material: dataset-specific narrative around the intestinal epithelium example.

## Why This Is A Separate Skill

- The notebook also contains differential cell-type abundance analysis, but that is a separate backend family with a different constructor, different required arguments, and different validation surface.
- DEG can be requested independently without any need to run DCT first.

## Notebook To Skill Mapping

- Dataset filtering to the control and Salmonella conditions became an input-contract example, not a hardcoded step.
- The Wilcoxon section became the default minimal execution pattern.
- The memento section became an alternate branch with explicit branch-specific kwargs and count-handling notes.
- Gene-specific violin plots were demoted to optional reporting because they do not define the reusable DEG contract.

## Example-Specific Material That Was Removed From The Trigger Surface

- The Haber dataset loader.
- The specific `TA` cell type.
- The example genes `Reg3b`, `Reg3g`, `Apoa1`, `Btf3`, `Serinc2`, and `Birc5`.
- The disease framing around one intestinal infection example.
