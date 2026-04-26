---
name: omicverse-single-cell-cellmatch-ontology
description: Map free-text cell-type annotations to the Cell Ontology (CL) via NLP-based sentence-transformer matching, with optional LLM-driven abbreviation expansion and Cell Taxonomy taxonomy enrichment. Use when standardising author-shorthand cell-type names to canonical CL terms (e.g. 'TIL-1' → 'tissue-resident memory CD8+ T cell, CL:0000625'), when running cross-cohort label harmonisation, or when reproducing `t_cellmatch`.
---

# OmicVerse Single-Cell — Cell Ontology Mapping (CellMatch)

## Goal

Take an annotated single-cell `AnnData` whose `obs[<cell_name_col>]` carries free-text cell-type labels (often author-specific shorthand) and **map every label to a canonical Cell Ontology (CL) term** via NLP. Output is `adata.obs['cell_ontology']` (matched CL term name), `adata.obs['cell_ontology_cl_id']` (e.g. `CL:0000084`), and `adata.obs['cell_ontology_score']` (cosine similarity).

Three escalating modes:
1. **Plain mapping** — sentence-transformer cosine-similarity match. Fastest; works on clean labels.
2. **LLM abbreviation expansion** — calls an LLM to expand `'TIL-1'` → `'tissue-resident memory CD8+ T cell'` *before* matching. Required when labels are short / acronym-heavy.
3. **Cell Taxonomy enhancement** — adds species + tissue context via the Jin 2023 Cell Taxonomy resource; matches against species-specific taxonomies as well as CL.

This skill addresses the gap noted in coverage analysis: `single-cell-annotation` skill handles CellTypist / SCSA / gpt4celltype annotators, but doesn't cover the *post-annotation* problem of mapping author-text to a standard ontology.

## Quick Workflow

### Mode 1 — Plain mapping (cleanest workflow when labels are full text)

1. Download the Cell Ontology to a local directory: `ov.single.download_cl(output_dir='new_ontology', filename='cl.json')`. ~10 MB; cached after first call.
2. Construct: `mapper = ov.single.CellOntologyMapper(cl_obo_file='new_ontology/cl.json', embeddings_path='new_ontology/ontology_embeddings.pkl', local_model_dir='./my_models')`. The first call computes ~3000-term embeddings (~30 s); subsequent calls load from `embeddings_path`.
3. **Map**: `mapping_results = mapper.map_adata(adata, cell_name_col='cell_label')`. Writes `obs['cell_ontology']`, `obs['cell_ontology_cl_id']`, `obs['cell_ontology_score']`.
4. **Summary**: `mapper.print_mapping_summary(mapping_results, top_n=15)` shows top-mapped labels and their match scores.
5. **Visualise**: `ov.pl.embedding(adata, basis='X_umap', color=['cell_label', 'cell_ontology'], wspace=0.55, ncols=2)` — side-by-side original and CL-mapped labels.

### Mode 2 — Plain mapping + LLM abbreviation expansion

6. Configure LLM access: `mapper.setup_llm_expansion(api_type='openai', api_key='sk-...', model='gpt-4o-2024-11-20', tissue_context='gut', species='mouse', study_context='Epithelial cells from small intestine and organoids of mice...')`. The contexts inform the LLM's expansion choices.
7. **Map with expansion**: `mapping_results = mapper.map_adata_with_expansion(adata=adata, cell_name_col='cell_label', threshold=0.5, expand_abbreviations=True)`. The LLM resolves `'TA-Early'` → `'transit amplifying cell, early'` first, then maps that to CL.
8. Use `api_type='custom_openai'` with `base_url=...` for OpenAI-compatible endpoints (Azure, Ollama, OhMyGPT, etc.).
9. **Cache**: abbreviation expansions are cached locally; subsequent runs reuse the cache (controlled by `cache_file` in `setup_llm_expansion`).

### Mode 3 — Cell Taxonomy enhancement (species + tissue aware)

10. Load the Cell Taxonomy resource (Jin 2023): `mapper.load_cell_taxonomy_resource('new_ontology/Cell_Taxonomy_resource.txt', species_filter=['Homo sapiens', 'Mus musculus'])`. ~30k entries; species filter avoids loading non-relevant taxa.
11. **Map with taxonomy**: `enhanced_results = mapper.map_adata_with_taxonomy(adata, cell_name_col='cell_label', new_col_name='enhanced_cell_ontology', expand_abbreviations=True, use_taxonomy=True, species='Mus musculus', tissue_context='Gut', threshold=0.3)`. Combines CL match + Cell Taxonomy match.
12. Visualise: `ov.pl.embedding(adata, color=['cell_label', 'cell_ontology', 'enhanced_cell_ontology', 'enhanced_cell_ontology_taxonomy_match', 'enhanced_cell_ontology_ct_id'], ...)` — full provenance.

## Interface Summary

```python
ov.single.download_cl(
    output_dir: str = 'new_ontology',
    filename: str = 'cl.json',
)
# Downloads CL to <output_dir>/<filename>; idempotent.

ov.single.CellOntologyMapper(
    cl_obo_file: str | None = None,
    embeddings_path: str | None = None,
    model_name: str = 'all-mpnet-base-v2',     # sentence-transformer
    local_model_dir: str | None = None,         # cache dir for the ST model
    auto_download: bool = True,
)
```

Core methods (all use cosine similarity at threshold `0.5` by default):
- `mapper.map_cells(cell_names: list, threshold: float = 0.5, use_llm_selection: bool = False, llm_candidates_count: int = 10) → dict` — for ad-hoc lookups (no AnnData).
- `mapper.map_adata(adata, cell_name_col=None, threshold=0.5, new_col_name='cell_ontology') → dict` — primary entry point; writes obs columns, returns the mapping result dict.
- `mapper.map_cells_with_expansion(cell_names, threshold, expand_abbreviations=True, tissue_context, species, study_context, use_llm_selection=True, llm_candidates_count=10)`.
- `mapper.map_adata_with_expansion(adata, cell_name_col=None, threshold=0.5, new_col_name='cell_ontology', expand_abbreviations=True, tissue_context=None, species=None, study_context=None, use_llm_selection=True, llm_candidates_count=10)`.
- `mapper.map_cells_with_taxonomy(cell_names, threshold, expand_abbreviations=True, use_taxonomy=True, species, tissue_context, study_context, use_llm_selection=True, llm_candidates_count=10)`.
- `mapper.map_adata_with_taxonomy(adata, cell_name_col=None, threshold=0.5, new_col_name='cell_ontology', expand_abbreviations=True, use_taxonomy=True, species=None, tissue_context=None, study_context=None)`.

Configuration & resources:
- `mapper.setup_llm_expansion(api_type='openai'|'custom_openai', api_key=None, model='gpt-3.5-turbo', base_url=None, cache_file='abbreviation_cache.json', tissue_context=None, species='human', study_context=None, extra_params=None)`.
- `mapper.load_cell_taxonomy_resource(taxonomy_file, species_filter=None)` — adds Cell Taxonomy as a second ontology.
- `mapper.load_ontology_mappings(popv_json_path)` — load PopV-style CL ID mappings.
- `mapper.create_ontology_resources(cl_obo_file, save_embeddings=True)` — re-build the embeddings index.
- `mapper.set_model(model_name, local_model_dir=None)` — switch sentence-transformer.
- `mapper.set_local_model(model_path)` — load from a local model dir (offline / air-gapped).
- `mapper.download_model()` — manual model download.

Inspection:
- `mapper.search_ontology_cells(keyword, case_sensitive=False, max_results=20)` — keyword search over CL terms.
- `mapper.list_ontology_cells(max_display=50, return_all=False)` — list CL terms.
- `mapper.browse_ontology_by_category(categories=None, max_per_category=10)` — browse by anatomy / function category.
- `mapper.find_similar_cells(cell_name, top_k=10)` — top-k CL matches to a single name.
- `mapper.find_similar_cells_taxonomy(cell_name, species=None, top_k=10)` — same but in Cell Taxonomy.
- `mapper.search_by_marker(markers, species=None, top_k=10)` — find CT terms by marker overlap.
- `mapper.get_cell_info(cell_name)` / `get_cell_info_taxonomy(cell_name, species=None)` — full CL / CT entry.
- `mapper.check_ontology_status()` — diagnostic: ontology size, embedding state, model state.
- `mapper.get_ontology_statistics()`, `mapper.get_statistics(mapping_results)`.
- `mapper.print_mapping_summary(mapping_results, top_n=10)`, `print_mapping_summary_taxonomy(...)`, `print_mapping_summary_with_ids(...)`.
- `mapper.show_expansion_summary(mapping_results)` — what the LLM expansion did.
- `mapper.test_abbreviation_detection(test_cases=None)` — sanity-check abbreviation detection.

Persistence:
- `mapper.save_embeddings(output_path=None)` / `load_embeddings(embeddings_path)`.
- `mapper.save_mapping_results(mapping_results, output_file)`.
- `mapper.clear_abbreviation_cache()`.

## Boundary

**Inside scope:**
- Plain CL mapping via cosine similarity.
- LLM-driven abbreviation expansion before matching.
- Cell Taxonomy resource enhancement for species + tissue context.
- Diagnostic helpers (search, list, browse, similar).
- Persistence of embeddings + abbreviation cache.

**Outside scope — separate skill:**
- *Building* the cell-type annotation in the first place — `single-cell-annotation` skill (CellTypist / SCSA / gpt4celltype).
- LLM-only annotation (no CL mapping) — `single-cell-annotation` skill's gptcelltype branch.
- Cross-modality label transfer (RNA → ATAC) — see `cross-modal-celltype-transfer` skill.
- Reference-based label transfer (atlas → query) — see `reference-label-transfer` skill.
- MetaTiME tumour-microenvironment cell-state annotation — see `single-cell-metatime-annotation`.

## Branch Selection

**Plain vs expansion vs taxonomy**
- **Plain (`map_adata`)**: labels are full descriptive text (e.g. `'CD8+ effector memory T cell'`). Don't pay the LLM cost.
- **Expansion (`map_adata_with_expansion`)**: labels include abbreviations (e.g. `'TA-Early'`, `'TIL-1'`, `'IEL'`). LLM expansion is critical for recall.
- **Taxonomy (`map_adata_with_taxonomy`)**: cohort is species-specific and tissue-specific; you want both CL and species-resolved Cell Taxonomy matches. Pass `species='Mus musculus'`, `tissue_context='Gut'`.

**`threshold` (cosine similarity cutoff)**
- 0.5 (default) — confidence floor; matches below this are reported as `'No clear match'`.
- 0.3 — looser; useful with abbreviation expansion where similarity scores are systematically lower.
- 0.7 — stricter; only high-confidence matches.

**`use_llm_selection` (in expansion mode)**
- `True` (default in `_with_expansion` / `_with_taxonomy`): LLM picks the best CL match from the top-`llm_candidates_count` cosine candidates. Higher accuracy, costs LLM tokens.
- `False`: pure cosine; no LLM at the matching step.

**`llm_candidates_count`**
- 10 (default) — gives the LLM enough alternatives to choose from.
- 5 — cheaper LLM calls; less robust.
- 20+ — overkill; the LLM rarely benefits from >15 candidates.

**`api_type` for `setup_llm_expansion`**
- `'openai'` — direct OpenAI API.
- `'custom_openai'` — OpenAI-compatible endpoint (Azure, Ollama, OhMyGPT, vLLM-OpenAI). Pass `base_url`.

**Sentence-transformer model**
- `'all-mpnet-base-v2'` (default) — best-quality general-purpose biomedical mapping at ~110M params.
- `'all-MiniLM-L6-v2'` — much smaller (~22M); 5–10× faster but ~5 % lower mapping accuracy.
- Switch with `mapper.set_model(model_name, local_model_dir)` or pass to constructor.

**Cache management**
- `embeddings_path='ontology_embeddings.pkl'` — first call writes; subsequent calls read.
- LLM expansion cache: `setup_llm_expansion(cache_file='abbreviation_cache.json')`. Persistent across sessions; clear with `mapper.clear_abbreviation_cache()`.
- Re-build embeddings when CL is updated: `mapper.create_ontology_resources(cl_obo_file, save_embeddings=True)`.

## Input Contract

- `AnnData` with `obs[cell_name_col]` populated; values are strings (cell-type labels). Empty / NaN entries skip mapping.
- For LLM expansion: outbound HTTPS access to OpenAI (or compatible); valid `api_key`. Costs ~$0.001 per cell-type label expanded; the cache makes repeated runs free.
- For Cell Taxonomy enhancement: download `Cell_Taxonomy_resource.txt` separately (URL in tutorial; ~50 MB) and pass to `load_cell_taxonomy_resource`.
- Sentence-transformer first run: ~110M-param model download (cached via HuggingFace; pass `local_model_dir` to control location).
- Hardware: CPU is sufficient (~30 s for ~3000-term embedding); GPU not required.

## Minimal Execution Patterns

```python
import omicverse as ov
import scanpy as sc

ov.plot_set()

# 1) Get the Cell Ontology
ov.single.download_cl(output_dir='new_ontology', filename='cl.json')

# 2) Construct mapper (first call computes embeddings; subsequent loads from cache)
mapper = ov.single.CellOntologyMapper(
    cl_obo_file='new_ontology/cl.json',
    embeddings_path='new_ontology/ontology_embeddings.pkl',
    local_model_dir='./my_models',
)

# 3) Plain mapping
mapping_results = mapper.map_adata(adata, cell_name_col='cell_label')
mapper.print_mapping_summary(mapping_results, top_n=15)

# 4) Inspect on UMAP
ov.pl.embedding(
    adata, basis='X_umap',
    color=['cell_label', 'cell_ontology', 'cell_ontology_cl_id'],
    wspace=0.55, ncols=2,
)
```

```python
# LLM-expansion mode (handles abbreviations)
mapper.setup_llm_expansion(
    api_type='openai', model='gpt-4o-2024-11-20',
    api_key='sk-...',
    tissue_context='gut',
    species='mouse',
    study_context='Epithelial cells from small intestine and organoids; some Salmonella or H. polygyrus infected.',
)

# OR via a custom OpenAI-compatible endpoint
# mapper.setup_llm_expansion(
#     api_type='custom_openai',
#     api_key='sk-...',
#     model='gpt-4.1-2025-04-14',
#     base_url='https://api.ohmygpt.com/v1',
# )

mapping_results = mapper.map_adata_with_expansion(
    adata=adata,
    cell_name_col='cell_label',
    threshold=0.5,
    expand_abbreviations=True,
)
mapper.print_mapping_summary(mapping_results, top_n=15)
mapper.show_expansion_summary(mapping_results)
```

```python
# Taxonomy-enhanced mode (species + tissue)
mapper.load_cell_taxonomy_resource(
    'new_ontology/Cell_Taxonomy_resource.txt',
    species_filter=['Homo sapiens', 'Mus musculus'],
)

enhanced_results = mapper.map_adata_with_taxonomy(
    adata,
    cell_name_col='cell_label',
    new_col_name='enhanced_cell_ontology',
    expand_abbreviations=True,
    use_taxonomy=True,
    species='Mus musculus',
    tissue_context='Gut',
    threshold=0.3,
)

ov.pl.embedding(
    adata, basis='X_umap',
    color=['cell_label', 'cell_ontology', 'enhanced_cell_ontology',
           'enhanced_cell_ontology_taxonomy_match',
           'enhanced_cell_ontology_ct_id'],
    wspace=0.55, ncols=2,
)
```

## Validation

- After `download_cl`: `cl.json` exists at the requested path; size ~10 MB.
- After construction: `mapper.check_ontology_status()` reports the ontology size (~3000 CL terms) and embedding state (loaded vs to-build).
- After `map_adata`: `obs['cell_ontology']` populated; `(obs['cell_ontology_score'] >= threshold).sum() / len(obs)` is the match rate. Healthy cohorts get >70 % at `threshold=0.5`. Lower rates suggest abbreviations — switch to expansion mode.
- After expansion mode: `mapper.show_expansion_summary(...)` tells you which labels got expanded; if zero labels were expanded, `setup_llm_expansion` wasn't called or the LLM rejected all.
- Cache freshness: re-running `map_adata_with_expansion` should hit the abbreviation cache for already-seen labels (no LLM cost). Watch the request count to confirm.
- Sanity check: pick a few labels manually and run `mapper.find_similar_cells(<label>, top_k=10)` — the top match should be biologically reasonable. If the top match is generic (e.g. `'native cell'`), the label is too vague.

## Resource Map

- See [`reference.md`](reference.md) for compact copy-paste snippets per mode.
- See [`references/source-grounding.md`](references/source-grounding.md) for verified `CellOntologyMapper` constructor + the ~30 method signatures, plus the docstring backfill log (the class-level docstring was a single emoji line; expanded to full).
- For the upstream cell-type annotation that produces the input `cell_label`, see existing skills `single-cell-annotation`, `reference-label-transfer`, `single-popv-annotation`, `omicverse-single-cell-metatime-annotation`.
- For multi-annotator consensus (after CL mapping), see `omicverse-single-cell-cellvote-consensus`.

## Examples
- "Map free-text cell labels in `obs['cell_label']` to Cell Ontology terms with the default mpnet-base-v2 sentence transformer."
- "Cohort uses heavy abbreviations like 'TIL-1' / 'TA-Early' — set up LLM expansion with `tissue_context='gut'` / `species='mouse'` and re-map."
- "Add Cell Taxonomy enhancement for a mouse-gut cohort with `threshold=0.3`."
- "List the top-10 CL matches for a single label without modifying any AnnData (`mapper.find_similar_cells('TIL-1', top_k=10)`)."

## References
- Tutorial notebook: [`t_cellmatch.ipynb`](https://omicverse.readthedocs.io/en/latest/Tutorials-single/t_cellmatch/) — full three-mode walkthrough on a mouse gut cohort.
- Cell Ontology: https://obofoundry.org/ontology/cl.html
- Cell Taxonomy resource: Jin *et al.* 2023.
- Live API verified — see [`references/source-grounding.md`](references/source-grounding.md).
