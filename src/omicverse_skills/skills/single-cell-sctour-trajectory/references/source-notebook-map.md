# Source Notebook Map

Source notebook:

- `omicverse_guide/docs/Tutorials-single/t_traj.ipynb`

## Partition

### Upstream Skills Or Shared Work

- preprocessing and PCA setup belong to the preprocessing skill
- `diffusion_map`, `slingshot`, `palantir`, PAGA, and Palantir follow-ups belong to the shared trajectory-inference skill

### This Skill

Notebook cells:

- restoring counts into `.X`
- `sc.pp.calculate_qc_metrics(...)`
- `TrajInfer(...).inference(method='sctour', alpha_recon_lec=0.5, alpha_recon_lode=0.5)`
- plotting `sctour_pseudotime`
- notebook-specific `1 - pseudotime` reversal

## Why This Is A Separate Skill

- the branch changes the backend family from graph-based trajectory inference to trainer-based latent dynamics
- the wrapper imports an external dependency not needed by the other trajectory branches
- it writes a different result family: `sctour_pseudotime`, `X_TNODE`, and `X_VF`

## Notebook-Specific Details Demoted Out Of The Main Trigger Surface

- pancreas-specific labels
- exact reconstruction-weight values
- the notebook's manual pseudotime reversal

Those remain examples, not the skill identity.
