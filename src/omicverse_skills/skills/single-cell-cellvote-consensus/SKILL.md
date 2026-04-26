---
name: omicverse-single-cell-cellvote-consensus
description: Multi-annotator consensus for single-cell labels via `ov.single.CellVote`. Combine labels from any subset of SCSA / gpt4celltype / GPTBioInsightor / scMulan / PopV per cluster; resolve disagreement either by LLM arbitration (online) or local-majority voting (offline). Output is `obs['CellVote_celltype']`. Use when you have multiple annotators on the same AnnData and need a single consensus label, or when reproducing `t_cellvote_pbmc3k`.
---

# OmicVerse Single-Cell — CellVote Multi-Annotator Consensus

## Goal

Take an annotated `AnnData` that already has labels from **two or more annotators** (e.g. `obs['scsa_annotation']`, `obs['gpt_celltype']`, `obs['gbi_celltype']`) and produce a per-cluster **consensus label** in `obs['CellVote_celltype']`. The class wraps the upstream annotators (so you can populate them via the same object) and arbitrates by either calling an LLM with the per-cluster candidate set + the cluster's marker genes (online; requires API key) or running a deterministic local-majority vote on lower-cased candidates (offline; tutorial-canonical demo).

CellVote is **OmicVerse's own** consensus layer; it's distinct from `single-popv-annotation` (POPV is a separate Bayesian voting tool from a different group). CellVote can include POPV as one of its inputs but is not a wrapper for POPV.

## Quick Workflow

1. Ensure each candidate annotator's labels are populated in `obs` columns (e.g. `scsa_annotation`, `gpt_celltype`, `gbi_celltype`). Either run the annotators yourself (see Branch Selection for the pre-baked `cv.scsa_anno()` / `cv.gpt_anno()` / `cv.gbi_anno()` etc.) or import labels from elsewhere.
2. Compute cluster-level marker genes — `marker_dict[cluster_id] = [top_genes...]`. CellVote uses these to give the LLM (or the local arbitrator) biological context per cluster.
3. **Construct**: `cv = CellVote(adata)`. Stores the AnnData reference; doesn't run anything yet.
4. **Vote** (online, LLM): `final_map = cv.vote(clusters_key='leiden', cluster_markers=marker_dict, celltype_keys=['scsa_annotation', 'gpt_celltype', 'gbi_celltype'], species='human', organization='PBMC', provider='openai', model='gpt-4o-mini', api_key='sk-...', result_key='CellVote_celltype')`. Calls the LLM once per cluster with the candidate labels + markers + species/organization context; LLM returns a single arbitrated label. Result lands in `adata.obs['CellVote_celltype']`.
5. **Vote (offline, local majority)**: monkey-patch `omicverse.single._cellvote.get_cluster_celltype = local_fn` *before* calling `cv.vote(...)`. The local function takes the same `(cluster_celltypes, cluster_markers, species, organization, model, base_url, provider, api_key)` signature and returns a dict mapping cluster → consensus label. The tutorial provides a one-line `pd.Series.value_counts().idxmax()` implementation that's deterministic and free.
6. **Inspect**: compare `obs[['leiden', 'scsa_annotation', 'gpt_celltype', 'gbi_celltype', 'CellVote_celltype']].head()` and the per-cluster summary `obs.groupby('leiden')[annot_cols].agg(lambda s: s.value_counts().index[0])`.

## Interface Summary

```python
ov.single.CellVote(adata: AnnData) -> None
# Constructor; stores adata reference.
```

Per-annotator runners (each populates a different `obs` column):
- `cv.scsa_anno()` — runs the SCSA pipeline; writes `obs['scsa_annotation']`. (Documented in `single-cell-annotation` skill.)
- `cv.gpt_anno()` — gpt4celltype-style LLM annotation; writes `obs['gpt_celltype']`.
- `cv.gbi_anno()` — GPTBioInsightor; writes `obs['gbi_celltype']`.
- `cv.scMulan_anno()` — scMulan large-language model; writes the corresponding column.
- `cv.popv_anno(ref_adata, ref_labels_key, ref_batch_key, query_batch_key=None, cl_obo_folder=None, save_path='tmp', prediction_mode='fast', methods=None, methods_kwargs=None)` — PopV consensus from a reference; writes its own column. (Use `single-popv-annotation` skill for stand-alone PopV.)

Consensus voting:
```python
cv.vote(
    clusters_key: str = None,
    cluster_markers: dict[str, list[str]] = None,
    celltype_keys: list[str] = [],
    model: str = 'gpt-3.5-turbo',
    base_url: str = None,
    species: str = 'human',
    organization: str = 'stomach',
    provider: str = 'openai',
    result_key: str = 'CellVote_celltype',
) -> dict[cluster_id, label]
# Writes adata.obs[result_key]; returns the cluster → label dict.
```

Offline arbitrator (monkey-patch pattern):
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
# Then cv.vote(...) uses the patched function.
```

## Boundary

**Inside scope:**
- Multi-annotator consensus arbitration (LLM or local-majority).
- The pre-baked annotator runners (`scsa_anno`, `gpt_anno`, `gbi_anno`, `scMulan_anno`, `popv_anno`).
- Custom arbitration via the monkey-patch hook.

**Outside scope — separate skill:**
- Running individual annotators stand-alone — see `single-cell-annotation` (CellTypist / SCSA / gpt4celltype) and `single-popv-annotation` (POPV).
- Cell Ontology mapping of the consensus labels — see `omicverse-single-cell-cellmatch-ontology`.
- MetaTiME tumour-microenvironment annotation — see `omicverse-single-cell-metatime-annotation`.
- Cross-modality label transfer — see `cross-modal-celltype-transfer`.

## Branch Selection

**Online vs offline arbitration**
- **Online (LLM)**: pass `provider='openai'` (or `'custom_openai'` with `base_url`) + `api_key`. The LLM sees the per-cluster candidate labels AND the marker genes AND the `species` / `organization` context — typically produces a more biologically informed consensus, especially when annotators disagree. Costs ~$0.001 per cluster on `gpt-4o-mini`; total run ~$0.01 for a 20-cluster cohort.
- **Offline (local majority)**: monkey-patch `cvmod.get_cluster_celltype` before voting. Deterministic, free, no network. Tutorial-recommended for CI / reproducibility.
- **Hybrid pattern**: run offline first to get a deterministic baseline; run online on a copy of the AnnData; compare. Disagreements highlight where biology is ambiguous.

**Annotator selection (`celltype_keys`)**
- **Two annotators**: still works — but ties are resolved alphabetically in local-majority mode (suboptimal). Add a third for tie-breaking.
- **Three annotators** (canonical): SCSA + gpt4celltype + GPTBioInsightor. Tutorial default.
- **More annotators**: monotonically improves consensus quality but increases LLM tokens and runtime.

**`provider` and `model`**
- `provider='openai'` + `model='gpt-4o-mini'` (tutorial demo): cheap, fast, accurate enough for consensus. Pass `api_key=...` from env.
- `provider='openai'` + `model='gpt-4o-2024-11-20'`: more expensive, more accurate on edge cases.
- `provider='custom_openai'` + `base_url='https://...'`: any OpenAI-compatible endpoint (Azure, Ollama, vLLM, OhMyGPT).

**`species` / `organization`**
- Required for the LLM context. `species='human'`, `organization='PBMC'` for the tutorial. Wrong values → LLM hallucinates non-existent cell types.
- For non-typical organisms, `organization` can be free-text (e.g. `'mouse spleen post-LPS'`) — the LLM uses it as prompt context.

**`cluster_markers` source**
- Pass cluster-level top markers (typically 5–10 per cluster). Build from `sc.tl.rank_genes_groups` results, or from your own DEG analysis.
- For LLM arbitration: more markers (5–10) help the LLM disambiguate; fewer (2–3) underdetermine the cluster.
- For local-majority: `cluster_markers` is *not* used (the arbitration is purely candidate-frequency); pass anyway since the API requires it.

**Annotator pre-population**
- If you haven't run annotators yet: use `cv.scsa_anno()` etc. to populate them via the same `cv` object.
- If labels come from elsewhere (a previous run, a collaborator's output): just write them into `obs[<keyname>]` manually before voting.
- Heterogeneous label vocabularies (e.g. one annotator says `'CD8 T'`, another says `'CD8+ T cell'`): **standardise first** with the cellmatch skill, then vote. Otherwise the local-majority arbitrator counts these as different votes.

## Input Contract

- `AnnData` with `obs[clusters_key]` populated (typically `'leiden'`).
- For each `key` in `celltype_keys`: `obs[key]` populated with cell-type strings. Empty / NaN values get treated as `'unknown'` candidate.
- `cluster_markers`: `dict[cluster_id, list[gene_str]]`. Cluster IDs must match the unique values in `obs[clusters_key]` (string types).
- For online voting: outbound HTTPS to OpenAI (or compatible); valid `api_key`. ~10s per cluster.
- For offline voting: monkey-patch in place before `cv.vote(...)`.

## Minimal Execution Patterns

```python
import os
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import omicverse as ov
from omicverse.single import CellVote

# 1) Load PBMC3k and preprocess
adata = sc.datasets.pbmc3k_processed()
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

# Build marker_dict for the LLM context
marker_dict = {
    cluster: list(adata.uns['rank_genes_groups']['names'][cluster][:5])
    for cluster in adata.obs['leiden'].cat.categories
}

# 2) Populate three annotator columns (here: simulated; in practice, run the annotators)
adata.obs['scsa_annotation'] = ...     # from cv.scsa_anno() or external
adata.obs['gpt_celltype']    = ...     # from cv.gpt_anno() or external
adata.obs['gbi_celltype']    = ...     # from cv.gbi_anno() or external
```

```python
# 3) Online consensus (requires API key)
cv = CellVote(adata)
final_map = cv.vote(
    clusters_key='leiden',
    cluster_markers=marker_dict,
    celltype_keys=['scsa_annotation', 'gpt_celltype', 'gbi_celltype'],
    species='human',
    organization='PBMC',
    provider='openai',
    model='gpt-4o-mini',
    # api_key='sk-...',
)
# adata.obs['CellVote_celltype'] populated; final_map is the cluster -> label dict
```

```python
# 4) Offline local-majority arbitration (deterministic, free)
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
    provider='openai', model='gpt-4o-mini',  # ignored by the patched function
)
print(adata.obs[['leiden', 'scsa_annotation', 'gpt_celltype',
                 'gbi_celltype', 'CellVote_celltype']].head())
```

```python
# 5) Per-cluster summary table
cols = ['leiden', 'scsa_annotation', 'gpt_celltype',
        'gbi_celltype', 'CellVote_celltype']
summary = (adata.obs
              .groupby('leiden')[cols[1:]]
              .agg(lambda s: s.value_counts().index[0]))
print(summary)
```

## Validation

- After `cv.vote(...)`: `adata.obs['CellVote_celltype']` populated with no NaN cells (clusters with no candidate labels become `'unknown'`).
- Compare `cv.vote(...)` keys against `obs[clusters_key].unique()` — every cluster should be in the returned dict.
- Disagreement diagnostic: for each cluster, count candidate labels — if all 3 annotators agree, consensus is trivially that label. Where they disagree, the LLM (online) typically picks the most marker-supported one; the local-majority (offline) just picks the most-frequent.
- LLM arbitrator failures: if `'unknown'` appears in many `CellVote_celltype` outputs, the LLM rejected the candidates (often because the cluster's marker genes don't match any of the candidate cell types). Inspect `cluster_markers[cluster]` and the candidate set manually.
- Heterogeneous vocabulary: if `'CD8 T'` and `'CD8+ T cell'` co-exist, local-majority sees them as different votes and may pick the third (different) annotator's answer. Standardise via the cellmatch skill before voting.
- For online runs: log the LLM cost per call (`gpt-4o-mini` is ~$0.001/cluster); total cost is `n_clusters × 0.001` USD.

## Resource Map

- See [`reference.md`](reference.md) for compact copy-paste snippets per mode.
- See [`references/source-grounding.md`](references/source-grounding.md) for verified `CellVote` constructor + method signatures and the LLM-arbitrator dispatch path.
- For individual annotators (CellTypist / SCSA / gpt4celltype), see existing `single-cell-annotation` skill.
- For POPV (a different consensus tool), see existing `single-popv-annotation` skill.
- For Cell Ontology mapping of CellVote outputs (label standardisation across cohorts), see `omicverse-single-cell-cellmatch-ontology`.

## Examples
- "Build a consensus label from `scsa_annotation` + `gpt_celltype` + `gbi_celltype` on PBMC3k via local-majority arbitration (offline)."
- "Run online CellVote with gpt-4o-mini on a 20-cluster cohort and report the cost."
- "Diagnose why CellVote keeps returning `'unknown'` for cluster 5 — likely the markers don't match any annotator's candidate."
- "Run all three annotators (`cv.scsa_anno()`, `cv.gpt_anno()`, `cv.gbi_anno()`), then `cv.vote(...)` for a one-shot consensus pipeline."

## References
- Tutorial notebook: [`t_cellvote_pbmc3k.ipynb`](https://omicverse.readthedocs.io/en/latest/Tutorials-single/t_cellvote_pbmc3k/) — PBMC3k offline-and-online consensus.
- Live API verified — see [`references/source-grounding.md`](references/source-grounding.md).
