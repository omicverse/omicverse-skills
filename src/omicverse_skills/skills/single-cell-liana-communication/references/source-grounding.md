# Source Grounding — LIANA cell-cell communication

## Interfaces Checked

`omicverse.single.run_liana`, `omicverse.single.to_comm_adata`, `omicverse.pl.ccc_heatmap`. Verified via `inspect.signature` + `inspect.getdoc` on each + direct reading of `omicverse/single/_liana.py` and `omicverse/pl/_ccc_heatmap.py`. Cross-checked against `t_ccc_liana.ipynb`.

## Live signatures

```python
ov.single.run_liana(
    adata, *,
    groupby: str,
    method: str = 'rank_aggregate',
    key_added: str = 'liana_res',
    inplace: bool = True,
    **kwargs                          # forwarded to liana.mt.<method>
)

ov.single.to_comm_adata(
    adata: anndata.AnnData | None = None, *,
    data: pd.DataFrame | None = None,
    result_uns_key: str | None = None,
    score_key: str = 'specificity_rank',
    pvalue_key: str = 'specificity_rank',
    inverse_score: bool = True,
    inverse_pvalue: bool = False,
    classification: str | Mapping[str, str] | None = None,
    classification_reference: str | pd.DataFrame | None = 'cellchat',
    classification_fallback: str | None = 'family',
    separator: str = '|',
) -> anndata.AnnData

ov.pl.ccc_heatmap(
    adata_or_comm, *,
    plot_type: str = 'dot',
    display_by: str = 'interaction',
    score_key: str = 'specificity_rank',
    pvalue_key: str = 'specificity_rank',
    classification_reference: str | pd.DataFrame = 'cellchat',
    classification_fallback: str = 'family',
    sender_use: str | list | None = None,
    receiver_use: str | list | None = None,
    signaling: str | list | None = None,
    pattern: str = 'incoming',
    pvalue_threshold: float = 0.05,
    top_n: int = 10,
    cmap: str = ...,
    figsize: tuple = ...,
    show: bool = False,
    ...
) -> (fig, ax)
```

## Source-grounded behavior

**`run_liana`:** thin wrapper that imports the `liana` PyPI package lazily (raises `ImportError` if missing); dispatches to `liana.mt.rank_aggregate` (or whatever `method` says) with the cohort's `adata` and `groupby` (cell-type column). All extra `**kwargs` flow through to LIANA's underlying call. Result is a long DataFrame keyed by (sender source, target receiver, ligand_complex, receptor_complex) plus the requested score columns. Convention: `inplace=True` writes to `adata.uns[key_added]`; `inplace=False` returns the DataFrame.

**`to_comm_adata`:** post-processes the long-format LIANA result into a directed-pair × interaction matrix:
- `obs_names` = `'sender|receiver'` (separator from `separator` kwarg).
- `var_names` = ligand-receptor complex IDs.
- `.X` is the score matrix (rows × interactions); `inverse_score=True` flips rank-style scores so larger = stronger.
- `.layers['pvalue']` carries the p-value matrix (with `inverse_pvalue` similarly).
- `var['classification']` is filled from the `classification_reference` (default CellChat pathway taxonomy); when an LR pair has no CellChat entry, falls back to `classification_fallback` (default `'family'` = chemokine / cytokine / TGF-β / etc.).
- Custom classification: pass a `pd.DataFrame` keyed by ligand_receptor or a `dict[lr_id → label]`.

**`ccc_heatmap`:** unified plotting entrypoint with eight `plot_type` modes; each mode resolves to a matplotlib figure. The `display_by` axis switches the underlying tensor between per-interaction and per-(sender, receiver) aggregated views. `sender_use` / `receiver_use` filter the directional axis; `signaling` filters by pathway label.

## Notebook ↔ skill alignment

| Notebook section | Skill section |
|---|---|
| `ov.pl.embedding(color='bulk_labels')` sanity check | Quick Workflow §2 |
| `ov.single.run_liana(method='rank_aggregate', resource_name='consensus')` | Quick Workflow §3 |
| `ov.single.to_comm_adata(classification_reference='cellchat', classification_fallback='family')` | Quick Workflow §4 |
| `ov.pl.ccc_heatmap` modes: dot / tile / heatmap (aggregation) / focused_heatmap / pathway_bubble / role_heatmap / role_network / sender_use / receiver_use | Quick Workflow §5-8; reference.md eight-mode block |
| Multi-condition stacked frame with `condition` column | Quick Workflow §9; reference.md multi-condition block |

## Docstring supplementation log

`run_liana` (21L), `to_comm_adata` (12L), `ccc_heatmap` (≥30L by inspection) — all already documented; left as-is for this skill.

## Reviewer-Run Empirical Checks

- All cited functions importable: `from omicverse.single import run_liana, to_comm_adata; from omicverse.pl import ccc_heatmap` ✓ (LIANA itself loaded lazily; `pip install liana` required).
- Tutorial cells use exactly the verified signatures — the `score_key='specificity_rank'` / `pvalue_key='specificity_rank'` symmetry is correct (LIANA reuses the same column for the rank-aggregate output).
- The eight `plot_type` modes were verified by reading the dispatch logic in `omicverse/pl/_ccc_heatmap.py`; each mode has its own dedicated rendering function.
- No live smoke run executed; the PBMC8k cohort + the bundled CellChat resource together would form a fast smoke target (~30 s for `run_liana`, sub-second for plotting).

## Caveats

- LIANA is an external dependency; first call after `pip install liana` may also pull resource files (cached locally after).
- The `rank_aggregate` method is a meta-method that runs CellPhoneDB / NATMI / Connectome / SingleCellSignalR / CellChat under the hood — runtime is the sum of those, typically minutes on a 10k-cell cohort.
- `to_comm_adata` returns a fresh AnnData; the original `adata` is not modified by this call (unlike `run_liana(inplace=True)` which writes to `uns`).
