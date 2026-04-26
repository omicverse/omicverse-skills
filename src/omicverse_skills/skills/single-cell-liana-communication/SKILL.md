---
name: omicverse-single-cell-liana-communication
description: LIANA+ ligand-receptor inference on single-cell AnnData via `ov.single.run_liana`, plus the OmicVerse cell-cell communication (CCC) plotting stack with `ov.pl.ccc_heatmap` and `ov.single.to_comm_adata`. Use when computing ligand-receptor scores from a labeled AnnData (`bulk_labels` / `cell_type` / etc.), when post-processing LIANA results into a CommAnnData for pathway-aware visualisation, or when reproducing `t_ccc_liana`.
---

# OmicVerse Single-Cell — LIANA+ Cell-Cell Communication

## Goal

Take a preprocessed annotated single-cell `AnnData` and run **LIANA+** ligand-receptor inference, producing per-(sender, receiver, ligand, receptor) score columns in `adata.uns['liana_res']`. Then post-process via `ov.single.to_comm_adata(...)` into a *communication AnnData* keyed on the sender→receiver direction, classified against the **CellChat / CellChatDB** pathway taxonomy. Visualise with `ov.pl.ccc_heatmap(...)` — eight `plot_type` modes covering dot / tile / aggregation heatmap / pathway-bubble / role-heatmap / role-network / focused-heatmap.

This skill is the LIANA-side complement to the existing `single-cell-cellphonedb-communication` skill. They produce comparable hit sets but LIANA aggregates across multiple base methods (CellPhoneDB, NATMI, Connectome, SingleCellSignalR, CellChat) into a `rank_aggregate` consensus — typically more robust than any single method.

## Quick Workflow

1. Load a preprocessed annotated `AnnData` (`obs[group_key]` populated; for the PBMC8k demo, `bulk_labels`).
2. Optional: visualise the cohort embedding to confirm cluster labels look right (`ov.pl.embedding(adata, basis='X_umap', color='bulk_labels', frameon='small')`).
3. **Run LIANA**: `ov.single.run_liana(adata, groupby='bulk_labels', method='rank_aggregate', resource_name='consensus', key_added='liana_res', inplace=True)`. Result lands at `adata.uns['liana_res']`.
4. **Pathway-aware reshape**: `comm_adata = ov.single.to_comm_adata(adata, result_uns_key='liana_res', score_key='specificity_rank', pvalue_key='specificity_rank', classification_reference='cellchat', classification_fallback='family')`. Returns a CommAnnData with one var per (ligand, receptor) and a `var['classification']` column mapping each pair to a CellChat pathway.
5. **Pathway dot plots**: `ov.pl.ccc_heatmap(adata, plot_type='dot', display_by='interaction', score_key='specificity_rank', ...)`. Multiple `plot_type` modes share a uniform interface — see Branch Selection below.
6. **Aggregation heatmap**: `plot_type='heatmap', display_by='aggregation'` shows sender→receiver totals.
7. **Sender / receiver focus**: pass `sender_use=<cluster>` or `receiver_use=<cluster>` to filter to a single direction.
8. **Pathway-focused / role-network views**: `plot_type='focused_heatmap'`, `'pathway_bubble'`, `'role_heatmap'`, `'role_network'` — pick by question (see Branch Selection).
9. **Multi-condition comparison**: stack LIANA results from multiple conditions, set `condition` column on the result frames, and re-run plotters with the multi-condition adata.

## Interface Summary

```python
ov.single.run_liana(
    adata, *,
    groupby: str,                           # required — the cell-type column
    method: str = 'rank_aggregate',          # or 'cellphonedb', 'natmi', 'connectome', 'sca'
    resource_name: str = 'consensus',        # ligand-receptor DB; 'consensus', 'cellphonedb', 'cellchat', etc.
    key_added: str = 'liana_res',
    inplace: bool = True,
    **kwargs                                 # forwarded to liana.mt.<method>
)
# adata.uns[key_added] -> pd.DataFrame with one row per (sender, receiver, ligand_complex, receptor_complex)
# columns include: 'specificity_rank', 'magnitude_rank', 'lr_means', 'cellphone_pvals', etc.
```

```python
ov.single.to_comm_adata(
    adata: AnnData | None = None, *,
    data: pd.DataFrame | None = None,
    result_uns_key: str | None = None,
    score_key: str = 'specificity_rank',
    pvalue_key: str = 'specificity_rank',
    inverse_score: bool = True,              # smaller rank = better → invert for visualisation
    inverse_pvalue: bool = False,
    classification: str | Mapping[str, str] | None = None,
    classification_reference: str | pd.DataFrame | None = 'cellchat',
    classification_fallback: str | None = 'family',
    separator: str = '|',
) -> AnnData
# Returns a "communication AnnData":
#   .obs   = (sender, receiver) directed pairs
#   .var   = (ligand, receptor) interactions; var['classification'] is the pathway label
#   .X     = score matrix (pairs × interactions)
#   .layers = {'pvalue': ...}
```

```python
ov.pl.ccc_heatmap(
    adata_or_comm, *,
    plot_type: str = 'dot',                  # see Branch Selection for full list
    display_by: str = 'interaction',         # 'interaction' or 'aggregation'
    score_key: str = 'specificity_rank',
    pvalue_key: str = 'specificity_rank',
    classification_reference='cellchat',
    classification_fallback='family',
    sender_use: str | list | None = None,
    receiver_use: str | list | None = None,
    signaling: str | list | None = None,     # focus on named pathways
    pattern: str = 'incoming',               # for role_heatmap / role_network: 'incoming' / 'outgoing'
    pvalue_threshold: float = 0.05,
    top_n: int = 10,
    cmap: str = ...,
    figsize: tuple = ...,
    show: bool = False,
    ...,
) -> (fig, ax)
```

## Boundary

**Inside scope:**
- LIANA fitting on a labeled `AnnData` with `run_liana`.
- Reshaping to `CommAnnData` with `to_comm_adata` (CellChat / CellPhoneDB / family classification).
- All eight `ov.pl.ccc_heatmap` modes.
- Sender / receiver / pathway focusing.
- Multi-condition comparison via stacked result frames with a `condition` column.

**Outside scope — separate skill:**
- CellPhoneDB-only inference + `CellChatViz` workflow — see `single-cell-cellphonedb-communication` (existing skill).
- Spatial cell-cell communication with COMMOT / FlowSig — see `spatial-tutorials` (existing skill).
- Building the upstream annotated AnnData — separate skill.
- Single-method LIANA (e.g. only NATMI) — supported by `method=` kwarg but rare in practice; the consensus `rank_aggregate` is the canonical choice.

## Branch Selection

**`method` for `run_liana`**
- `'rank_aggregate'` (default) — consensus across CellPhoneDB / NATMI / Connectome / SingleCellSignalR / CellChat. **The right default**.
- `'cellphonedb'` — single-method permutation-based; matches the CellPhoneDB v3 paper.
- `'natmi'` — specificity edge-weighting; less conservative than CellPhoneDB.
- `'connectome'` — scaled mean expression; gives both `magnitude` and `specificity` ranks.
- `'sca'` (SingleCellSignalR) — LRscore-based; conservative.

**`resource_name`**
- `'consensus'` (default) — union of CellChat + CellPhoneDB + Ramilowski + others.
- `'cellphonedb'` — strict CellPhoneDB v4 only.
- `'cellchat'` — strict CellChat only.
- For organism-specific runs, append `_human` / `_mouse` to the resource name.

**`score_key` and `pvalue_key`**
- `'specificity_rank'` — both score and pvalue (LIANA reuses the same column for the rank-aggregate output). For other methods, use the method-specific columns (`'cellphone_pvals'` for CellPhoneDB, `'lr_logfc'` for NATMI, etc.).
- `'magnitude_rank'` — alternate score that emphasises raw expression magnitude over specificity.
- `inverse_score=True` (default in `to_comm_adata`): converts rank-style scores (lower=better) to score-style (higher=better) for visualisation.

**`plot_type` for `ccc_heatmap`**
- `'dot'` — sender × receiver bubble plot, dot size = score, colour = pvalue. Default exploratory view.
- `'tile'` — heat tile: rows = interactions, cols = sender→receiver pairs.
- `'heatmap'` (with `display_by='aggregation'`) — sender × receiver aggregated total CCC strength.
- `'focused_heatmap'` — restrict to a named pathway (`signaling=['ECM/Adhesion']`).
- `'pathway_bubble'` — pathway × sender→receiver bubble; needs `signaling=[...]`.
- `'role_heatmap'` — sender / receiver role intensity per cluster (`pattern='incoming' | 'outgoing'`).
- `'role_network'` — directed network plot of dominant signalling roles.

**`display_by`**
- `'interaction'` — per (ligand, receptor) detail.
- `'aggregation'` — summed across interactions per (sender, receiver). Use with `plot_type='heatmap'`.

**`sender_use` / `receiver_use`**
- `sender_use='CD34+'` — show only signals originating from CD34+. Useful for "where does CD34+ talk to?"
- `receiver_use='CD34+'` — show only signals received by CD34+. Useful for "what regulates CD34+?"

**`classification_reference`**
- `'cellchat'` (default) — CellChat pathway taxonomy.
- `'family'` (fallback) — chemokine / cytokine / TGF-β / etc. families.
- Pass a custom `pd.DataFrame` keyed by ligand-receptor for arbitrary pathway labels.

**Multi-condition comparison**: build a stacked frame with a `condition` column, then call `ov.pl.ccc_heatmap` with that frame. Allows side-by-side baseline vs. stimulated visualisation.

## Input Contract

- `AnnData` with `obs[group_key]` populated (categorical or strings); typical group_keys: `'bulk_labels'`, `'cell_type'`, `'leiden'`.
- `adata.X` log-normalised expression. LIANA expects this scale; passing raw counts gives nonsensical specificity scores.
- `pip install liana` (the omicverse install does not ship LIANA by default).
- For `to_comm_adata` with `classification_reference='cellchat'`: requires the CellChat ligand-receptor mapping shipped under `omicverse.external.cellchat` (always available in OmicVerse).

## Minimal Execution Patterns

```python
import omicverse as ov
import scanpy as sc
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

ov.plot_set(font_path='Arial')

# 1) Load a labeled AnnData
adata = sc.datasets.pbmc8k()                # or your own
adata.obs['bulk_labels'] = ...               # cluster labels

ov.pl.embedding(adata, basis='X_umap', color='bulk_labels', frameon='small')

# 2) Run LIANA
ov.single.run_liana(
    adata,
    groupby='bulk_labels',
    method='rank_aggregate',
    resource_name='consensus',
    key_added='liana_res',
    inplace=True,
)
print(adata.uns['liana_res'].head())

# 3) Reshape into a CommAnnData with CellChat pathway classification
comm_adata = ov.single.to_comm_adata(
    adata,
    result_uns_key='liana_res',
    score_key='specificity_rank',
    pvalue_key='specificity_rank',
    classification_reference='cellchat',
    classification_fallback='family',
)
print(comm_adata)
print(comm_adata.var['classification'].value_counts().head())

# 4) Dot plot — default exploratory view
fig, ax = ov.pl.ccc_heatmap(
    adata,
    plot_type='dot',
    display_by='interaction',
    score_key='specificity_rank',
    pvalue_key='specificity_rank',
    classification_reference='cellchat',
    classification_fallback='family',
    show=False,
)
```

```python
# Sender / receiver focused dot plots
fig, ax = ov.pl.ccc_heatmap(
    comm_adata, plot_type='dot', display_by='interaction',
    score_key='specificity_rank', pvalue_key='specificity_rank',
    sender_use='CD34+', top_n=6, pvalue_threshold=0.05,
    show=False,
)
fig, ax = ov.pl.ccc_heatmap(
    adata, plot_type='dot', display_by='interaction',
    score_key='specificity_rank', pvalue_key='specificity_rank',
    receiver_use='CD34+', top_n=5, pvalue_threshold=0.05,
    show=False,
)

# Aggregation heatmap (sender × receiver totals)
fig, ax = ov.pl.ccc_heatmap(
    adata, plot_type='heatmap', display_by='aggregation',
    score_key='specificity_rank', pvalue_key='specificity_rank',
    classification_reference='cellchat', classification_fallback='family',
    cmap='YlGnBu', figsize=(4, 3), show=False,
)

# Focus on a specific pathway
focus = comm_adata.var['classification'].dropna().astype(str).unique().tolist()
focus = [v for v in focus if v not in {'Unclassified', 'nan'}]
focus_pathway = focus[0] if focus else 'Unclassified'

fig, ax = ov.pl.ccc_heatmap(
    adata, plot_type='focused_heatmap',
    signaling=[focus_pathway],
    min_interaction_threshold=0.0,
    cmap='YlGnBu', figsize=(4, 3), show=False,
)

# Pathway-bubble for a specific pathway
fig, ax = ov.pl.ccc_heatmap(
    adata, plot_type='pathway_bubble',
    signaling=['ECM/Adhesion'],
    top_n=5, figsize=(3, 5), show=False,
)

# Role heatmap (incoming vs outgoing)
fig, ax = ov.pl.ccc_heatmap(
    adata, plot_type='role_heatmap',
    pattern='incoming', cmap='Greens',
    figsize=(4, 3), show=False,
)

# Role network for a pathway
fig, ax = ov.pl.ccc_heatmap(
    adata, plot_type='role_network',
    signaling=['ECM/Adhesion'],
    cmap='Greens', figsize=(8, 5), show=False,
)
```

## Validation

- After `run_liana`: `adata.uns['liana_res']` is a non-empty DataFrame; key columns include `source` (sender), `target` (receiver), `ligand_complex`, `receptor_complex`, `specificity_rank`, `magnitude_rank`. For PBMC8k, expect 100k+ rows (n_clusters² × n_LR_pairs).
- After `to_comm_adata`: `comm_adata.shape` is roughly `(n_directional_pairs, n_LR_pairs)`. `comm_adata.var['classification']` should have most entries non-null; null entries are unmapped LR pairs.
- For `plot_type='dot'`: dot size and colour scale should be readable; if the figure is empty, the `pvalue_threshold` is too strict OR `score_key` was wrong.
- Sender/receiver focus: `sender_use=<cluster>` should reduce the number of bubbles substantially; if the figure is identical to the unfocused one, the cluster name was misspelled.
- `plot_type='focused_heatmap'` empty: the named pathway has no LR pairs above the threshold OR the pathway name doesn't exist (check `comm_adata.var['classification'].unique()`).
- For multi-condition plots: each condition's frame must have the same column schema (`source`, `target`, `ligand_complex`, `receptor_complex`, `specificity_rank` minimally).
- Compare LIANA `rank_aggregate` hits against CellPhoneDB-only hits on the same cohort: typically 60–80 % overlap at the top-100; large disagreement (<30 %) means one method has bad calibration on this cohort.

## Resource Map

- See [`reference.md`](reference.md) for compact copy-paste snippets.
- See [`references/source-grounding.md`](references/source-grounding.md) for verified `run_liana` / `to_comm_adata` / `ccc_heatmap` signatures.
- For the alternative CellPhoneDB-only workflow + CellChatViz, see `single-cell-cellphonedb-communication`.
- For spatial cell-cell communication (COMMOT / FlowSig), see `spatial-tutorials`.

## Examples
- "Run LIANA `rank_aggregate` on a PBMC8k cohort grouped by `bulk_labels` and produce a sender × receiver dot plot."
- "Reshape the LIANA result into a CommAnnData classified against CellChat pathways and plot the role network for ECM/Adhesion."
- "Focus the dot plot on signals originating from CD34+ at `pvalue_threshold=0.05`."
- "Compare LIANA results between baseline and stimulated conditions on the same cohort."

## References
- Tutorial notebook: [`t_ccc_liana.ipynb`](https://omicverse.readthedocs.io/en/latest/Tutorials-single/t_ccc_liana/) — PBMC8k LIANA + CellChat-style visualisation.
- LIANA+ paper: Dimitrov *et al.* 2022, *Nature Communications* — "Comparison of methods and resources for cell-cell communication inference".
- Live API verified — see [`references/source-grounding.md`](references/source-grounding.md).
