# LIANA cell-cell communication — quick commands

```python
import omicverse as ov
import scanpy as sc
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
ov.plot_set(font_path='Arial')

# Labeled AnnData (replace with yours)
adata = sc.datasets.pbmc8k()
# adata.obs['bulk_labels'] should be populated

# 1) Run LIANA — rank_aggregate consensus across multiple methods
ov.single.run_liana(
    adata,
    groupby='bulk_labels',
    method='rank_aggregate',
    resource_name='consensus',
    key_added='liana_res',
    inplace=True,
)
print(adata.uns['liana_res'].head())
```

## Reshape into pathway-classified CommAnnData

```python
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
```

## Eight `ccc_heatmap` modes

```python
import matplotlib.pyplot as plt

# 1. Default dot plot — sender × receiver bubbles
fig, ax = ov.pl.ccc_heatmap(
    adata, plot_type='dot', display_by='interaction',
    score_key='specificity_rank', pvalue_key='specificity_rank',
    classification_reference='cellchat', classification_fallback='family',
    show=False,
)

# 2. Tile heatmap (interactions × pairs)
fig, ax = ov.pl.ccc_heatmap(
    adata, plot_type='tile', display_by='interaction',
    score_key='specificity_rank', pvalue_key='specificity_rank',
    top_n=10, pvalue_threshold=0.05,
    cmap='viridis', figsize=..., show=False,
)

# 3. Aggregation heatmap (sender × receiver totals)
fig, ax = ov.pl.ccc_heatmap(
    adata, plot_type='heatmap', display_by='aggregation',
    score_key='specificity_rank', pvalue_key='specificity_rank',
    classification_reference='cellchat', classification_fallback='family',
    cmap='YlGnBu', figsize=(4, 3), show=False,
)

# 4. Pathway-focused heatmap
focus_pathway = (
    comm_adata.var['classification'].dropna().astype(str)
    .replace('Unclassified', None).dropna().iloc[0]
)
fig, ax = ov.pl.ccc_heatmap(
    adata, plot_type='focused_heatmap',
    signaling=[focus_pathway],
    min_interaction_threshold=0.0,
    cmap='YlGnBu', figsize=(4, 3), show=False,
)

# 5. Pathway bubble
fig, ax = ov.pl.ccc_heatmap(
    adata, plot_type='pathway_bubble',
    signaling=['ECM/Adhesion'],
    top_n=5, figsize=(3, 5), show=False,
)

# 6. Role heatmap — incoming or outgoing
fig, ax = ov.pl.ccc_heatmap(
    adata, plot_type='role_heatmap',
    pattern='incoming',                 # or 'outgoing'
    cmap='Greens', figsize=(4, 3), show=False,
)

# 7. Role network for a specific pathway
fig, ax = ov.pl.ccc_heatmap(
    adata, plot_type='role_network',
    signaling=['ECM/Adhesion'],
    cmap='Greens', figsize=(8, 5), show=False,
)

# 8. Sender / receiver focus
fig, ax = ov.pl.ccc_heatmap(
    comm_adata, plot_type='dot', display_by='interaction',
    score_key='specificity_rank', pvalue_key='specificity_rank',
    sender_use='CD34+', top_n=6, pvalue_threshold=0.05, show=False,
)
fig, ax = ov.pl.ccc_heatmap(
    adata, plot_type='dot', display_by='interaction',
    score_key='specificity_rank', pvalue_key='specificity_rank',
    receiver_use='CD34+', top_n=5, pvalue_threshold=0.05, show=False,
)
```

## Multi-condition stack

```python
import numpy as np

# Stack two conditions with a 'condition' column
frames = []
for condition, delta in [('baseline', 0.0), ('stimulated', 0.03)]:
    frame = adata.uns['liana_res'].copy()
    frame['condition'] = condition
    frame['specificity_rank'] = np.clip(frame['specificity_rank'] + delta, 0, 1)
    frames.append(frame)
adata_by_context = adata.copy()
adata_by_context.uns['liana_res'] = pd.concat(frames, ignore_index=True)

# Pass through the same plot_type='dot' / 'heatmap' machinery
fig, ax = ov.pl.ccc_heatmap(
    adata_by_context, plot_type='dot', display_by='interaction',
    score_key='specificity_rank', pvalue_key='specificity_rank',
    show=False,
)
```
