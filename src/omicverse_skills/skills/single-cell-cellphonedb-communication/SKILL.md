---
name: omicverse-single-cell-cellphonedb-communication
description: Analyze single-cell cell-cell communication with OmicVerse CellPhoneDB and CellChat-style visualization. Use when converting an OmicVerse CellPhoneDB notebook into a reusable skill, when running CellPhoneDB ligand-receptor analysis on annotated AnnData, or when choosing pathway aggregation, layout, signaling-role, and bubble-plot branches for downstream communication summaries.
---

# OmicVerse Single-Cell CellPhoneDB Communication

## Goal

Run the reusable CellPhoneDB communication spine on single-cell `AnnData`: verify the expression matrix and cell-type annotations, run `ov.single.run_cellphonedb_v5(...)`, convert the results into the visualization-ready interaction `AnnData`, then optionally continue into aggregated network plots, pathway summaries, ligand-receptor contribution views, bubble plots, chord diagrams, and signaling-role analysis through `ov.pl.CellChatViz(...)`. Keep the skill centered on one processed CellPhoneDB result object rather than on the tutorial dataset.

## Quick Workflow

1. Inspect the input `AnnData`, especially the cell-type column, gene-symbol convention, and whether the matrix looks compatible with CellPhoneDB.
2. Ensure the CellPhoneDB database archive is available or let the wrapper download it automatically.
3. Run `ov.single.run_cellphonedb_v5(...)` with explicit filtering, permutation, and output settings.
4. Treat the returned `cpdb_results` dict and `adata_cpdb` object as the shared handoff point for all downstream visualization branches.
5. Initialize `viz = ov.pl.CellChatViz(adata_cpdb, palette=...)`.
6. Choose the downstream branch the user actually asked for: aggregated network, pathway-level aggregation, ligand-receptor extraction, bubble/chord views, or signaling-role analysis.
7. Validate the expected layers, sender/receiver annotations, and pathway metadata before trusting any plot.

## Interface Summary

- `ov.single.run_cellphonedb_v5(adata, cpdb_file_path, celltype_key='celltype', min_cell_fraction=0.005, min_genes=200, min_cells=3, iterations=1000, threshold=0.1, pvalue=0.05, threads=10, output_dir=None, temp_dir=None, cleanup_temp=True, debug=False, separator='|', **kwargs)` runs CellPhoneDB statistical analysis and returns both the raw result dict and a visualization-ready interaction `AnnData`.
- The wrapper filters cell types by `min_cell_fraction`, then runs `scanpy` filtering with `min_genes` and `min_cells` before invoking CellPhoneDB.
- `ov.pl.CellChatViz(adata_cpdb, palette=None)` initializes the CellChat-style visualization object for the processed communication result.
- `compute_aggregated_network(pvalue_threshold=0.05, use_means=True)` returns interaction-count and interaction-weight matrices for sender-receiver cell types.
- `compute_pathway_communication(method='mean', min_lr_pairs=1, min_expression=0.1)` is a real branch point with aggregation options `mean`, `sum`, `max`, and `median`.
- `get_significant_pathways_v2(pathway_communication=None, strength_threshold=0.1, pvalue_threshold=0.05, min_significant_pairs=1)` summarizes and ranks pathway-level communication.
- `netVisual_aggregate(signaling, layout='circle', ...)` branches on `layout='circle'` or `layout='hierarchy'`.
- `netVisual_chord_cell(signaling=None, group_celltype=None, sources=None, targets=None, ..., normalize_to_sender=True)` branches across pathway selection, sender/receiver filtering, optional grouping, and arc normalization.
- `extractEnrichedLR(signaling, pvalue_threshold=0.05, mean_threshold=0.1, min_cell_pairs=1, geneLR_return=False)` returns pathway-contributing ligand-receptor pairs.
- `netVisual_bubble_marsilea(..., signaling=None, sources_use=None, targets_use=None, group_pathways=True, transpose=False, scale=None, ...)` branches on pathway grouping, sender/receiver filters, transposition, and scaling mode.
- `netAnalysis_computeCentrality(signaling=None, pvalue_threshold=0.05, use_weight=True)` computes sender/receiver/mediator-style centrality summaries.
- `netAnalysis_signalingRole_heatmap(pattern='outgoing', signaling=None, row_scale=True, ..., min_threshold=0.1)` branches on `pattern='outgoing'` or `pattern='incoming'`.

## Boundary

- Keep this as one skill because the notebook contains one tight communication-analysis job with multiple downstream views on the same `adata_cpdb` result.
- Do not split aggregated-network, pathway, bubble, chord, and signaling-role analysis into separate skills by default, because they all require the same processed CellPhoneDB interaction object and are weak on their own without it.
- Split only when the user already has a valid CellPhoneDB visualization `AnnData` and explicitly wants a visualization-only helper without rerunning analysis.
- Do not absorb unrelated preprocessing, annotation, or trajectory workflows here.

## Branch Selection

- Use `min_cell_fraction` to exclude rare cell types before CellPhoneDB runs; increase it when noisy small groups are diluting the interaction graph.
- Use `iterations` high enough for a valid CellPhoneDB statistical run. Very small values can fail in the underlying CellPhoneDB implementation; a bounded smoke path should still keep permutations in a valid range.
- Use `cleanup_temp=True` for the normal path; switch it off only when you need to inspect temporary input files.
- Use `debug=True` only when you intentionally want CellPhoneDB intermediate tables.
- Use `method='mean'` in `compute_pathway_communication` for the notebook-style pathway summary; switch to `sum`, `max`, or `median` only when you intentionally want a different pathway aggregation rule.
- Use `layout='circle'` in `netVisual_aggregate` for the notebook-style summary view; use `layout='hierarchy'` only when sender/receiver partitioning is the main question.
- Use `signaling=None` to summarize all pathways, or a string/list of pathways when you want pathway-specific views.
- Use `sources` and `targets` in chord views, or `sources_use` and `targets_use` in bubble views, when the user wants a sender- or receiver-focused slice.
- Use `group_celltype` in chord views only when several cell types should collapse into a higher-level group.
- Use `scale='row'`, `scale='column'`, `scale='row_minmax'`, or `scale='column_minmax'` in the bubble view only when the user explicitly wants normalized visual comparison instead of raw communication strength.
- Use `pattern='outgoing'` to study sender roles and `pattern='incoming'` to study receiver roles in signaling-role heatmaps.

## Input Contract

- Start from an annotated `AnnData` with a valid cell-type column.
- Prefer human gene symbols for CellPhoneDB analysis, because CellPhoneDB is built around human ligand-receptor resources.
- Ensure the expression matrix is compatible with CellPhoneDB; the tutorial path uses log-normalized values and explicitly avoids scaled data.
- Expect the visualization-ready output object to store `means` and `pvalues` as layers and `sender` and `receiver` as observation metadata.
- Expect pathway-aware downstream methods to rely on interaction classification metadata in the variable table.

## Minimal Execution Patterns

```python
import omicverse as ov

cpdb_results, adata_cpdb = ov.single.run_cellphonedb_v5(
    adata,
    cpdb_file_path=cpdb_db,
    celltype_key="cell_type",
    iterations=1000,
    pvalue=0.05,
)

viz = ov.pl.CellChatViz(adata_cpdb, palette=palette)
count_matrix, weight_matrix = viz.compute_aggregated_network(
    pvalue_threshold=0.05,
    use_means=True,
)
```

```python
pathway_comm = viz.compute_pathway_communication(
    method="mean",
    min_lr_pairs=2,
    min_expression=0.1,
)
sig_pathways, pathway_summary = viz.get_significant_pathways_v2(
    pathway_comm,
    strength_threshold=0.1,
    pvalue_threshold=0.05,
    min_significant_pairs=1,
)
fig, ax = viz.netVisual_aggregate(
    signaling=sig_pathways[:1],
    layout="circle",
)
```

```python
enriched_lr = viz.extractEnrichedLR(
    signaling=["Signaling by Fibroblast growth factor"],
    pvalue_threshold=0.05,
    mean_threshold=0.1,
)
h = viz.netVisual_bubble_marsilea(
    sources_use=None,
    targets_use=None,
    signaling=["Signaling by Fibroblast growth factor"],
    scale=None,
)
centrality_scores = viz.netAnalysis_computeCentrality(
    signaling=None,
    pvalue_threshold=0.05,
    use_weight=True,
)
```

## Validation

- Check that the cell-type column exists before calling `ov.single.run_cellphonedb_v5(...)`.
- Check that the returned `cpdb_results` includes at least `means` and `pvalues`.
- Check that `adata_cpdb.layers` contains `means` and `pvalues`.
- Check that `adata_cpdb.obs` contains `sender` and `receiver`.
- Check that pathway-aware analyses have interaction classification metadata in `adata_cpdb.var`.
- If `compute_pathway_communication(...)` returns no pathways, inspect whether filtering thresholds removed too many interactions.
- If the user requests pathway-specific plots, confirm that the requested pathway names actually exist before plotting.
- If you only validated a bounded local smoke path, say so explicitly and do not claim a full publication-scale reproduction.

## Resource Map

- Read the branch-selection reference when choosing CellPhoneDB filtering, pathway aggregation, bubble scaling, or signaling-role branches.
- Read the source-grounding reference before extending the skill with more interface-specific claims.
- Read the notebook-map reference when deciding whether a future communication notebook belongs here or should become a downstream-only visualization skill.
- Read the compatibility reference when CellPhoneDB permutation settings, pathway significance thresholds, or optional plotting dependencies matter.
