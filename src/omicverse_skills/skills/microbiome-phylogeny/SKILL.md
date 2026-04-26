---
name: omicverse-microbiome-phylogeny
description: Build an ASV-level phylogenetic tree (MAFFT alignment + FastTree GTR+Γ), attach it to a 16S AnnData, and run phylogenetically-aware diversity (Faith PD, weighted / unweighted UniFrac). Use when adding tree-aware metrics on top of the basic 16S amplicon AnnData, or when the cohort needs UniFrac-driven beta-diversity instead of Bray-Curtis.
---

# OmicVerse Microbiome — Phylogenetic Diversity

## Goal

Take a 16S amplicon `AnnData` (samples × ASVs, with ASV centroid sequences in `var`) plus the ASV centroid FASTA from the amplicon pipeline, build a per-cohort phylogenetic tree (MAFFT alignment → FastTree under GTR+Γ), attach it to `adata.uns['tree']`, and run the **phylogenetically-aware** diversity stack: **Faith's PD** (alpha) and **unweighted / weighted UniFrac** (beta). Compare against non-phylogenetic baselines (Bray-Curtis) to confirm that phylogeny adds signal — or that it doesn't.

This skill assumes the upstream `omicverse-microbiome-16s-amplicon-dada2` skill has already produced the AnnData and the ASV FASTA (the FASTA is the `nochim['asv']` output of the stepwise pipeline, or is regenerable from the saved AnnData).

## Quick Workflow

1. Load the saved AnnData; identify the ASV FASTA from the workdir.
2. **Build the tree**: `tree = ov.alignment.build_phylogeny(asvs_fasta, workdir, mafft_threads=8, fasttree_threads=4)`. Returns a dict with `aligned`, `tree`, `newick` paths/strings.
3. **Attach** to AnnData: `ov.micro.attach_tree(adata, newick=tree['newick'])`. Stores at `adata.uns['tree']`; tip-name validation runs by default.
4. **Faith PD (alpha)**: include `'faith_pd'` in `Alpha(adata).run(metrics=[...])`. Result column lands in `adata.obs['faith_pd']`.
5. **UniFrac (beta)**: `Beta(adata).run(metric='unweighted_unifrac')` and `metric='weighted_unifrac'`. Distance matrices land in `adata.obsp[<metric>]`. Both call out to the `unifrac` package which expects rooted newick (FastTree midpoint-roots automatically).
6. **Compare with Bray-Curtis**: run all three on the same AnnData; PCoA on each; visually compare PC1/PC2 cluster separation. UniFrac usually agrees with Bray-Curtis on cohorts dominated by abundance shifts; it diverges when the *which* taxa shifted matters more than how much.
7. Save the now-tree-augmented AnnData (`adata.write_h5ad(...)`).

## Interface Summary

Tree construction:
```python
ov.alignment.build_phylogeny(
    asvs_fasta: str,
    workdir: Optional[str] = None, *,
    mafft_mode: str = 'auto',         # MAFFT alignment strategy (auto / FFT-NS-2 / L-INS-i / E-INS-i)
    fasttree_model: str = 'gtr',      # 'gtr' (default; nucleotide GTR) or 'jc'
    gamma: bool = True,               # GTR+Γ rate variation
    mafft_threads: int = 4,
    fasttree_threads: Optional[int] = None,
    overwrite: bool = False,
) -> Dict[str, str]
# returns {'aligned': <path>, 'tree': <path to .nwk>, 'newick': <tree text>, ...}
```

Attaching the tree:
```python
ov.micro.attach_tree(
    adata: AnnData,
    newick: Optional[str] = None,
    tree_path: Optional[Union[str, Path]] = None,
    prune: bool = True,
    store_key: str = 'tree',
    strict: bool = False,
) -> AnnData
```

Phylogenetically-aware metrics (existing `ov.micro` API; require the tree to be attached):
- Alpha: `Alpha(adata).run(metrics=['faith_pd', ...])` — computes Faith's phylogenetic diversity (sum of branch lengths on the subtree induced by present features) per sample. Reads tree from `adata.uns[tree_key]` (default `'tree'`).
- Beta: `Beta(adata).run(metric='unweighted_unifrac' | 'weighted_unifrac', rarefy=...)`. UniFrac distances read the same `adata.uns[tree_key]`.
- Both fall back to a clear `ValueError` if the tree is missing.

## Boundary

**Inside scope:**
- ASV-level tree construction with MAFFT + FastTree.
- Tree attachment + tip-name validation against `adata.var_names`.
- Faith PD alpha, unweighted/weighted UniFrac beta.
- Comparison plots / metrics across phylo and non-phylo metrics.

**Outside scope — separate skill:**
- Reference-based tree placement (e.g., placing query ASVs on a SILVA backbone via `gappa` / `epang`) — current implementation builds a de-novo tree from the cohort.
- IQ-TREE / RAxML replacement of FastTree — `build_phylogeny` is FastTree-only; IQ-TREE would be a separate function.
- Cross-cohort tree harmonisation — see `omicverse-microbiome-meta-analysis` for the meta-analysis-side caveats; this skill builds one tree per cohort.
- Differential abundance with tree-aware tests (gneiss, balance trees) — not implemented in `ov.micro`.
- 18S / ITS phylogeny — implementation assumes 16S rRNA (uses GTR+Γ; OK for 16S, less ideal for fast-evolving fungal regions).

## Branch Selection

**`mafft_mode`**
- `'auto'` (default) — MAFFT picks based on sequence count: `FFT-NS-2` for >200 sequences, `L-INS-i` otherwise. Right answer 95 % of the time.
- `'L-INS-i'` — most accurate for <500 sequences; slower. Use when the cohort is small and you care about alignment quality.
- `'E-INS-i'` — for sequences with very long unalignable regions; rare in 16S V3-V4 amplicons.

**`fasttree_model='gtr'` + `gamma=True`** is the canonical 16S choice. Switching to JC69 (`'jc'`) is for didactic purposes; never report a real result on JC.

**`fasttree_threads`** — `None` (FastTree decides) is fine. Increasing gives diminishing returns past ~4 threads on typical 16S cohorts.

**Faith PD is depth-sensitive** — like `observed_otus`, it grows with sequencing depth. Always rarefy first OR pair with Shannon (which is depth-corrected in expectation).

**Unweighted vs weighted UniFrac**
- Unweighted: presence/absence; sensitive to rare taxa, sensitive to differences in *which* taxa are present. Use when the question is "do these communities have different membership?"
- Weighted: abundance-weighted; less sensitive to rare taxa, more sensitive to differences in *how much* of each taxon is present. Use when the question is "how different are the communities in their dominant members?"
- Report both — they often disagree, and the disagreement is informative.
- Don't rarefy weighted UniFrac (`rarefy=False` default for UniFrac in `ov.micro`); rarefy unweighted UniFrac if depths are very uneven.

**`attach_tree(prune=True)`** is the right default: prunes tree tips not in `adata.var_names`, which speeds up UniFrac. `strict=True` raises on tip mismatches; useful in CI.

## Input Contract

- `adata` from the 16S amplicon skill: `var_names` are ASV ids (`asv00001`, etc.), `var['asv_seq']` (or similar) carries the centroid sequence string.
- `asvs_fasta`: path to the centroid FASTA (`nochim['asv']` from the stepwise pipeline or directly from `workdir/asv/asvs.fasta` after `amplicon_16s_pipeline`). Tip names must equal `var_names`.
- External binaries on `$PATH`: `mafft`, `FastTreeMP`/`FastTree`. Both are installable via `conda install -c bioconda mafft fasttree`.
- For UniFrac: `pip install unifrac` (or `conda install -c bioconda unifrac`). The `Beta` wrapper imports lazily and raises a clear error if missing.
- `adata.obs[group_col]` populated for any group-coloured PCoA plots (not required for tree construction itself).

## Minimal Execution Patterns

```python
import omicverse as ov
import anndata as ad
from pathlib import Path

ov.plot_set()

BASE = Path('/scratch/.../cache/16s/run_mothur_sop')
adata = ad.read_h5ad(BASE / 'mothur_sop_16s.h5ad')

# 1) Build the tree from the ASV FASTA
tree = ov.alignment.build_phylogeny(
    asvs_fasta=str(BASE / 'asv' / 'asvs.fasta'),
    workdir=str(BASE / 'phylogeny'),
    mafft_threads=8,
    fasttree_threads=4,
)
print('tree file:', tree['tree'])

# 2) Attach (prunes tips not in adata.var_names by default)
ov.micro.attach_tree(adata, newick=tree['newick'])
print("tree tips:", len(adata.uns['tree']))

# 3) Faith PD alpha (alongside Shannon / observed)
ov.micro.Alpha(adata).run(metrics=['shannon', 'observed_otus', 'faith_pd'])
print(adata.obs[['shannon', 'observed_otus', 'faith_pd']].describe())

# 4) UniFrac beta
b = ov.micro.Beta(adata)
b.run(metric='unweighted_unifrac', rarefy=False)
b.run(metric='weighted_unifrac',   rarefy=False)
b.run(metric='braycurtis', rarefy=True)            # baseline for comparison

# 5) PCoA on each metric
import matplotlib.pyplot as plt, pandas as pd
fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
for ax, dist_key, title in zip(
    axes,
    ['braycurtis', 'unweighted_unifrac', 'weighted_unifrac'],
    ['Bray-Curtis', 'unweighted UniFrac', 'weighted UniFrac'],
):
    ord_ = ov.micro.Ordinate(adata, dist_key=dist_key); ord_.pcoa(n=2)
    pct = ord_.proportion_explained() * 100.0
    coords = pd.DataFrame(adata.obsm[f'{dist_key}_pcoa'],
                          index=adata.obs_names, columns=['PC1', 'PC2'])
    for g, sub in adata.obs.groupby('group'):
        ax.scatter(coords.loc[sub.index, 'PC1'], coords.loc[sub.index, 'PC2'],
                   label=g, s=70, alpha=0.85, edgecolor='k')
    ax.set_xlabel(f'PC1 ({pct[0]:.1f}%)')
    ax.set_ylabel(f'PC2 ({pct[1]:.1f}%)')
    ax.set_title(title)
axes[-1].legend(title='group', frameon=True)
plt.tight_layout(); plt.show()

# 6) Save the tree-augmented AnnData
adata.write_h5ad(BASE / 'mothur_sop_16s_tree.h5ad')
```

## Validation

- `tree['tree']` exists on disk (`.nwk`); `tree['newick']` is non-empty text starting with `(`.
- After `attach_tree`: `adata.uns['tree']` is a non-empty newick string; `adata.uns['micro']['tree_tips']` matches expected ASV count (default `prune=True` so tips are at most `adata.n_vars`).
- After Faith PD: `adata.obs['faith_pd']` has no NaNs and grows monotonically with `observed_otus` on rarefied data (a sanity check).
- UniFrac matrix shapes: `adata.obsp['unweighted_unifrac'].shape == (n_obs, n_obs)`; off-diagonal mean is in [0, 1] (UniFrac distance is bounded).
- The two UniFrac variants typically *disagree* on which samples are most similar — that disagreement is itself the signal.
- If `build_phylogeny` errors: `mafft` / `FastTree` not on `$PATH`. Check `which mafft` / `which FastTreeMP` first.
- If `Beta(...).run(metric='unweighted_unifrac')` errors: `unifrac` package not installed (`pip install unifrac`) OR tip names don't match `adata.var_names` (rerun `attach_tree(strict=True)` to surface mismatches).

## Resource Map

- See [`reference.md`](reference.md) for compact copy-paste snippets.
- See [`references/source-grounding.md`](references/source-grounding.md) for verified `build_phylogeny` / `attach_tree` / UniFrac dispatch and the docstring backfill log.
- For the AnnData ingest that produces the input, see `omicverse-microbiome-16s-amplicon-dada2`.
- For DA on tree-augmented data (Wilcoxon / pyDESeq2 / ANCOM-BC), see `omicverse-microbiome-da-comparison` — none of them use the tree directly, but reporting alongside Faith PD / UniFrac is common.

## Examples
- "Build a FastTree GTR+Γ phylogeny from `asvs.fasta`, attach to the AnnData, and run Faith PD + unweighted/weighted UniFrac."
- "Compare Bray-Curtis, unweighted UniFrac, and weighted UniFrac PCoAs side by side on the mothur SOP cohort."
- "Confirm UniFrac tip names match `adata.var_names` (`attach_tree(strict=True)`)."
- "Save the tree-augmented AnnData for downstream meta-analysis."

## References
- Tutorial notebook: [`t_16s_phylogeny.ipynb`](https://omicverse.readthedocs.io/en/latest/Tutorials-microbiome/t_16s_phylogeny/) — full tree → Faith PD → UniFrac walkthrough on the mothur SOP.
- Live API verified — see [`references/source-grounding.md`](references/source-grounding.md).
