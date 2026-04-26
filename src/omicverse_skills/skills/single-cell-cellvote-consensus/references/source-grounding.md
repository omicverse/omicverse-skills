# Source Grounding — CellVote consensus

## Interfaces Checked

`omicverse.single.CellVote` and the module-level `omicverse.single._cellvote.get_cluster_celltype` arbitration hook. Verified via `inspect.signature` + `inspect.getdoc` + direct reading of `omicverse/single/_cellvote.py`. Cross-checked against `t_cellvote_pbmc3k.ipynb`.

## Live signatures

```python
ov.single.CellVote(adata) -> None

# Per-annotator runners (each writes a fixed obs column):
.scsa_anno(self)
.gpt_anno(self)
.gbi_anno(self)
.scMulan_anno(self)
.popv_anno(self, ref_adata, ref_labels_key, ref_batch_key,
           query_batch_key=None, cl_obo_folder=None,
           save_path='tmp', prediction_mode='fast',
           methods=None, methods_kwargs=None)

# Consensus arbitration:
.vote(self,
      clusters_key=None,
      cluster_markers=None,
      celltype_keys=[],
      model='gpt-3.5-turbo',
      base_url=None,
      species='human',
      organization='stomach',
      provider='openai',
      result_key='CellVote_celltype',
)
# Returns dict[cluster_id, label]; writes adata.obs[result_key].
```

Module-level hook (monkey-patchable for offline mode):
```python
omicverse.single._cellvote.get_cluster_celltype(
    cluster_celltypes: dict[cluster_id, list[str]],
    cluster_markers: dict[cluster_id, list[str]],
    species: str,
    organization: str,
    model: str,
    base_url: str | None,
    provider: str,
    api_key: str | None = None,
    **kwargs,
) -> dict[cluster_id, label]
```

## Source-grounded behavior

**Constructor:** stores `self.adata = adata` (reference, not copy). All methods mutate `self.adata.obs` in place.

**Per-annotator runners** delegate to the corresponding annotator module:
- `scsa_anno` → `pySCSA` → writes `obs['scsa_annotation']`.
- `gpt_anno` → `gpt4celltype` → writes `obs['gpt_celltype']`.
- `gbi_anno` → `GPTBioInsightor` → writes `obs['gbi_celltype']`.
- `scMulan_anno` → scMulan model → writes the corresponding column.
- `popv_anno` → POPV with the supplied reference; writes the consensus column.

Each runner is a thin convenience: stand-alone usage of these annotators is documented in `single-cell-annotation` and `single-popv-annotation`.

**`vote(...)` flow:**
1. Build `cluster_celltypes: dict[cluster_id, list[label]]` by collecting `obs[key]` for each `key` in `celltype_keys`, grouped by `obs[clusters_key]`. Within a cluster, takes the modal label from each annotator (so each annotator contributes one vote per cluster).
2. Calls `omicverse.single._cellvote.get_cluster_celltype(cluster_celltypes, cluster_markers, species, organization, model, base_url, provider, api_key)`. The default implementation is an LLM call.
3. The function returns `dict[cluster_id, final_label]`.
4. The result is broadcast into `obs[result_key]` per cell (all cells in cluster X get the cluster X consensus label).

**Monkey-patch pattern:** because `vote(...)` looks up `get_cluster_celltype` at module level (not as a method), replacing it before calling `vote(...)` redirects all subsequent calls to the custom function. After the offline run, `importlib.reload(cvmod)` restores the LLM-based original.

**LLM arbitrator** (default `get_cluster_celltype`):
- Builds a prompt with: per-cluster candidate labels, top-5 marker genes, species, organization context.
- Calls the named provider/model via the chat-completions API.
- Parses the LLM's structured output into `dict[cluster_id, label]`.
- Errors / parse failures fall back to `'unknown'` for the affected cluster.
- `provider='openai'` uses the OpenAI Python client; `provider='custom_openai'` requires `base_url=...` and routes through the same client.

## Notebook ↔ skill alignment

| Notebook section | Skill section |
|---|---|
| Synthesise three noisy annotator columns from canonical markers | Quick Workflow §1-2 (mirrored as "populate annotators") |
| Monkey-patch `cvmod.get_cluster_celltype = local_majority_arbitration` | Quick Workflow §5; Branch Selection (offline) |
| `cv.vote(clusters_key='leiden', cluster_markers=marker_dict, celltype_keys=[...], species='human', organization='PBMC', provider='openai', model='gpt-4o-mini')` | Quick Workflow §4 (online) and §5 (offline; same call after patch) |
| Per-cluster summary `obs.groupby('leiden')[cols[1:]].agg(...)` | Quick Workflow §6; reference.md summary block |
| `RUN_ONLINE = False` toggle (online only when API key available) | Branch Selection (online vs offline) |

## Docstring supplementation log

`CellVote` class (6L) — adequate; documents the multi-backend nature.
Methods (`vote` 29L, `scsa_anno` 9L, `gpt_anno` 10L, `gbi_anno` 13L, `popv_anno` 23L, `scMulan_anno` 10L) — all documented Numpy-style.

No supplementation done in this skill's pass.

## Reviewer-Run Empirical Checks

- All cited functions importable: `from omicverse.single import CellVote; import omicverse.single._cellvote as cvmod` ✓
- Monkey-patch pattern verified by reading the source: `vote()` calls `cvmod.get_cluster_celltype(...)` (module-level lookup), so replacing the module attribute before invocation is the correct pattern.
- Default `model='gpt-3.5-turbo'` and `organization='stomach'` are placeholders — tutorial overrides both (`model='gpt-4o-mini'`, `organization='PBMC'`).
- Offline arbitrator (the tutorial-canonical local-majority implementation) takes the full `(cluster_celltypes, cluster_markers, species, organization, model, base_url, provider, api_key)` signature even though it only reads `cluster_celltypes` — the module-level dispatch passes all kwargs through.
- No live smoke run executed; the offline mode is a pure-Python deterministic test that any CI environment can run.
