# Source Notebook Map

Source notebook:

- `omicverse_guide/docs/Tutorials-single/t_traj.ipynb`

## Partition

### Reused Upstream Skill

Notebook cells:

- dataset load
- `ov.pp.preprocess(...)`
- `ov.pp.scale(...)`
- `ov.pp.pca(...)`
- PCA variance plotting
- initial UMAP display

Reason:

- this is the same stable preprocessing job already covered by the preprocessing skill

### This Skill: Shared `TrajInfer` Family

Notebook cells:

- `TrajInfer(...)`
- `set_origin_cells(...)`
- `set_terminal_cells(...)`
- `inference(method='diffusion_map')`
- `inference(method='slingshot', ...)`
- `inference(method='palantir', ...)`
- `ov.utils.cal_paga(...)`
- `ov.utils.plot_paga(...)`
- `palantir_cal_branch(...)`
- `palantir_cal_gene_trends(...)`
- `palantir_plot_gene_trends(...)`

Reason:

- one shared input contract
- one shared wrapper object
- realistic users often ask for these branches as alternative trajectory methods on the same prepared object

### Separate Skill: `sctour`

Notebook cells:

- raw-count reset into `.X`
- `sc.pp.calculate_qc_metrics(...)`
- `TrajInfer(...).inference(method='sctour', ...)`
- `adata.obs['sctour_pseudotime'] = 1 - ...`

Reason:

- different backend family
- trainer-based execution
- extra dependency
- different compute profile

## Notebook-Specific Details Demoted Out Of The Main Trigger Surface

- pancreas development dataset choice
- `Ductal`, `Alpha`, `Beta`, `Delta`, `Epsilon` labels
- marker genes such as `Pax4` and `Ins2`
- exact notebook figure titles

Those stay as examples, not as the skill identity.
