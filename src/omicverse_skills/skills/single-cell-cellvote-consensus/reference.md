# CellVote consensus — quick commands

## Setup — annotator columns + cluster markers

```python
import pandas as pd
import scanpy as sc
import anndata as ad
import omicverse as ov
from omicverse.single import CellVote

# Build top-N marker_dict for the LLM (or local) arbitrator
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
marker_dict = {
    cluster: list(adata.uns['rank_genes_groups']['names'][cluster][:5])
    for cluster in adata.obs['leiden'].cat.categories
}

# Each annotator must be present in obs[<key>]
# Either run them via the cv methods OR import labels from elsewhere
```

## Pre-baked annotator runners (CellVote-managed)

```python
cv = CellVote(adata)
cv.scsa_anno()                # writes obs['scsa_annotation']
cv.gpt_anno()                 # writes obs['gpt_celltype']
cv.gbi_anno()                 # writes obs['gbi_celltype']
# cv.scMulan_anno()           # optional
# cv.popv_anno(ref_adata, ref_labels_key, ref_batch_key, ...)   # optional
```

## Online consensus (LLM arbitrator)

```python
cv = CellVote(adata)
final_map = cv.vote(
    clusters_key='leiden',
    cluster_markers=marker_dict,
    celltype_keys=['scsa_annotation', 'gpt_celltype', 'gbi_celltype'],
    species='human', organization='PBMC',
    provider='openai',          # or 'custom_openai' with base_url
    model='gpt-4o-mini',         # or 'gpt-4o-2024-11-20'
    # api_key='sk-...',           # if not in env
    result_key='CellVote_celltype',
)
print(adata.obs[['leiden', 'CellVote_celltype']].head())
print(final_map)                # dict: cluster -> label
```

## Offline local-majority arbitration (deterministic, free)

```python
import omicverse.single._cellvote as cvmod

def local_majority_arbitration(cluster_celltypes, cluster_markers,
                               species, organization, model, base_url,
                               provider, api_key=None, **kwargs):
    out = {}
    for cl, cand in cluster_celltypes.items():
        if not cand:
            out[cl] = 'unknown'
        else:
            s = pd.Series(cand).str.lower()
            out[cl] = s.value_counts().idxmax()
    return out

cvmod.get_cluster_celltype = local_majority_arbitration

cv = CellVote(adata)
final_map_offline = cv.vote(
    clusters_key='leiden',
    cluster_markers=marker_dict,
    celltype_keys=['scsa_annotation', 'gpt_celltype', 'gbi_celltype'],
    species='human', organization='PBMC',
    provider='openai', model='gpt-4o-mini',  # ignored by the patched fn
)
```

## Per-cluster summary

```python
cols = ['leiden', 'scsa_annotation', 'gpt_celltype',
        'gbi_celltype', 'CellVote_celltype']
print(adata.obs[cols].head())

summary = (adata.obs
              .groupby('leiden')[cols[1:]]
              .agg(lambda s: s.value_counts().index[0]))
print(summary)
```

## Disagreement diagnostic

```python
disagreement_mask = (
    (adata.obs['scsa_annotation'].astype(str)
       != adata.obs['gpt_celltype'].astype(str))
    | (adata.obs['gpt_celltype'].astype(str)
       != adata.obs['gbi_celltype'].astype(str))
)
print(f'{disagreement_mask.sum()} / {len(adata)} cells have disagreeing annotators')
print(adata.obs.loc[disagreement_mask,
                    ['leiden', 'scsa_annotation', 'gpt_celltype',
                     'gbi_celltype', 'CellVote_celltype']].head(20))
```

## Online + offline comparison (side-by-side)

```python
# Restore the original arbitrator after the offline run
import importlib
importlib.reload(cvmod)

cv2 = CellVote(adata.copy())
result_online = cv2.vote(
    clusters_key='leiden',
    cluster_markers=marker_dict,
    celltype_keys=['scsa_annotation', 'gpt_celltype', 'gbi_celltype'],
    species='human', organization='PBMC',
    provider='openai', model='gpt-4o-mini',
)

# Compare
import pandas as pd
compare = pd.DataFrame({
    'offline': pd.Series(final_map_offline),
    'online':  pd.Series(result_online),
})
print(compare[compare['offline'] != compare['online']])
```
