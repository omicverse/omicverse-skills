# Source Grounding — Phylogenetic Diversity

## Interfaces Checked

`omicverse.alignment.build_phylogeny` and `omicverse.micro.attach_tree`. Verified via `inspect.signature` + `inspect.getdoc` + direct reading of `omicverse/alignment/build_phylogeny.py`, `omicverse/micro/_phylo.py`, `omicverse/micro/_diversity.py` (Alpha + Beta UniFrac dispatch). Cell-by-cell cross-checking against `t_16s_phylogeny.ipynb`.

## Live signatures

```python
ov.alignment.build_phylogeny(
    asvs_fasta: str,
    workdir: Optional[str] = None, *,
    mafft_mode: str = 'auto',
    fasttree_model: str = 'gtr',
    gamma: bool = True,
    mafft_threads: int = 4,
    fasttree_threads: Optional[int] = None,
    overwrite: bool = False,
) -> Dict[str, str]

ov.micro.attach_tree(
    adata: AnnData,
    newick: Optional[str] = None,
    tree_path: Optional[Union[str, Path]] = None,
    prune: bool = True,
    store_key: str = 'tree',
    strict: bool = False,
) -> AnnData
```

## Source-grounded behavior

**`build_phylogeny`:**
- Runs MAFFT first to align the ASV centroids; output FASTA at `workdir/aligned.fasta`. `mafft_mode='auto'` lets MAFFT pick `FFT-NS-2` for large input, `L-INS-i` otherwise.
- Calls FastTree (preferring `FastTreeMP` when the binary is available) on the alignment with `-nt -gtr -gamma` flags by default (matching `fasttree_model='gtr', gamma=True`).
- Strips vsearch-style `;size=N` annotations from the tip labels in the resulting newick so that `unifrac` and `ete3` can parse them directly.
- Returns `{'aligned': <fasta path>, 'tree': <nwk path>, 'newick': <tree text>, 'workdir': <workdir>}`.
- Errors gracefully when `mafft` / `FastTree` not on `$PATH` — surfaces the missing binary name in the message.

**`attach_tree`:**
- Stores the newick text at `adata.uns[store_key]` (default `'tree'`) and the tip count at `adata.uns['micro']['tree_tips']`.
- `prune=True` (default): drops tree tips not in `adata.var_names` before storage; speeds up downstream UniFrac substantially on cohorts where ASVs have been filtered post-tree-build.
- `strict=True`: raises `ValueError` if any tree tip is not in `adata.var_names`. Useful in CI; default `False` is a soft warning so the function can be used on partially-pruned cohorts.
- Pass `newick` (string) or `tree_path` (file path) — exactly one.

**Faith PD (`Alpha.run(metrics=['faith_pd'])`):**
- Reads tree from `adata.uns[tree_key]` (default `'tree'`).
- Implementation uses `skbio.diversity.alpha.faith_pd`; Faith PD = sum of branch lengths in the subtree induced by present taxa per sample.
- Depth-sensitive — see Quick Workflow / Branch Selection notes; rarefy first or report alongside Shannon.

**UniFrac (`Beta.run(metric='unweighted_unifrac' | 'weighted_unifrac')`):**
- Lazily imports the `unifrac` package; raises a clear `ImportError` if missing.
- Reads tree from `adata.uns[tree_key]`. Tip names must match `adata.var_names` exactly (this is why `attach_tree(prune=True)` is the right default).
- Distance matrix returned as a `pd.DataFrame` indexed by `adata.obs_names`, also written to `adata.obsp[<metric>]` when `write_to_obsp=True` (default).
- `rarefy=None` resolves to `True` only for `unweighted_unifrac` and Bray-Curtis / Jaccard; `False` for `weighted_unifrac` (since weighted is intended to be depth-stable).

## Notebook ↔ skill alignment

| Notebook section | Skill section |
|---|---|
| Load saved `mothur_sop_16s.h5ad` and locate `asvs.fasta` | Quick Workflow §1 |
| `build_phylogeny(mafft_threads=8, fasttree_threads=4)` | Quick Workflow §2; Branch Selection (mafft_mode, fasttree_model) |
| `attach_tree(adata, newick=tree['newick'])` + `tree_tips` print | Quick Workflow §3; Validation |
| `Alpha(adata).run(metrics=['shannon', 'observed_otus', 'faith_pd'])` + 3-panel boxplot | Quick Workflow §4; reference.md alpha block |
| `Beta(adata).run(metric='unweighted_unifrac' | 'weighted_unifrac' | 'braycurtis')` + off-diag magnitude print | Quick Workflow §5-6; reference.md three-panel PCoA |
| Save `mothur_sop_16s_tree.h5ad` | Quick Workflow §7 |

## Docstring supplementation log

`build_phylogeny` (concise but complete; documents return shape and the `;size=N` stripping behavior), `attach_tree` (18L Numpy-style with `prune` / `strict` semantics), `Alpha.run` (5L), `Beta.run` (17L) — all already documented; no supplementation needed.

`Alpha.shannon` / `Alpha.observed` / `Beta.braycurtis` / `Ordinate.nmds` / `Ordinate.proportion_explained` were filled as part of the `ov.micro` docstring audit (see 16S amplicon skill's `source-grounding.md`); those don't directly affect this skill's claims.

## Reviewer-Run Empirical Checks

- All cited functions importable: `from omicverse.alignment import build_phylogeny; from omicverse.micro import attach_tree, Alpha, Beta, Ordinate` ✓
- Tutorial cells use exactly the verified signatures.
- `;size=N` stripping behaviour confirmed in the source — newick output is `unifrac`-compatible without further processing.
- No live smoke run executed; the mothur SOP demo cohort is the intended smoke target (requires `mafft` + `FastTree` + `unifrac` on the test environment, which the project's CI sets up via conda).
