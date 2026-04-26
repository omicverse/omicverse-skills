# Source Grounding — Cell Ontology mapping (CellMatch)

## Interfaces Checked

`omicverse.single.CellOntologyMapper` and `omicverse.single.download_cl`. Verified via `inspect.signature` + `inspect.getdoc` + direct reading of `omicverse/single/_cellmatch.py`. Cross-checked against `t_cellmatch.ipynb`.

## Live signatures

```python
ov.single.download_cl(
    output_dir: str = 'new_ontology',
    filename: str = 'cl.json',
)

ov.single.CellOntologyMapper(
    cl_obo_file: str | None = None,
    embeddings_path: str | None = None,
    model_name: str = 'all-mpnet-base-v2',
    local_model_dir: str | None = None,
    auto_download: bool = True,
)
```

## Public methods (28 listed below; signatures verified by `inspect`)

```python
# Configuration
.setup_llm_expansion(api_type='openai', api_key=None, model='gpt-3.5-turbo',
                     base_url=None, cache_file='abbreviation_cache.json',
                     tissue_context=None, species='human',
                     study_context=None, extra_params=None)
.set_model(model_name, local_model_dir=None)
.set_local_model(model_path)
.download_model()

# Resources
.load_ontology_mappings(popv_json_path)
.load_cell_taxonomy_resource(taxonomy_file, species_filter=None)
.create_ontology_resources(cl_obo_file, save_embeddings=True)

# Inspection
.search_ontology_cells(keyword, case_sensitive=False, max_results=20)
.list_ontology_cells(max_display=50, return_all=False)
.browse_ontology_by_category(categories=None, max_per_category=10)
.find_similar_cells(cell_name, top_k=10)
.find_similar_cells_taxonomy(cell_name, species=None, top_k=10)
.search_by_marker(markers, species=None, top_k=10)
.get_cell_info(cell_name)
.get_cell_info_taxonomy(cell_name, species=None)
.get_ontology_statistics()
.check_ontology_status()
.test_abbreviation_detection(test_cases=None)

# Core mapping (returns dict; map_adata variants also write obs columns)
.map_cells(cell_names, threshold=0.5, use_llm_selection=False, llm_candidates_count=10)
.map_cells_with_expansion(cell_names, threshold=0.5, expand_abbreviations=True,
                          tissue_context=None, species=None, study_context=None,
                          use_llm_selection=True, llm_candidates_count=10)
.map_cells_with_taxonomy(cell_names, threshold=0.5, expand_abbreviations=True,
                         use_taxonomy=True, species=None, tissue_context=None,
                         study_context=None, use_llm_selection=True,
                         llm_candidates_count=10)
.map_adata(adata, cell_name_col=None, threshold=0.5, new_col_name='cell_ontology')
.map_adata_with_expansion(adata, cell_name_col=None, threshold=0.5,
                          new_col_name='cell_ontology', expand_abbreviations=True,
                          tissue_context=None, species=None, study_context=None,
                          use_llm_selection=True, llm_candidates_count=10)
.map_adata_with_taxonomy(adata, cell_name_col=None, threshold=0.5,
                         new_col_name='cell_ontology', expand_abbreviations=True,
                         use_taxonomy=True, species=None, tissue_context=None,
                         study_context=None)

# Reporting
.expand_abbreviations(cell_names, force_expand=False, save_cache=True,
                      tissue_context=None, species=None, study_context=None)
.print_mapping_summary(mapping_results, top_n=10)
.print_mapping_summary_taxonomy(mapping_results, top_n=10)
.print_mapping_summary_with_ids(mapping_results, top_n=10)
.show_expansion_summary(mapping_results)
.get_statistics(mapping_results)

# Persistence
.save_embeddings(output_path=None)
.load_embeddings(embeddings_path)
.save_mapping_results(mapping_results, output_file)
.clear_abbreviation_cache()
```

## Source-grounded behavior

**Constructor + `download_cl`:**
- `download_cl(output_dir, filename)` writes `cl.json` from the OBO Foundry mirror; idempotent (skipped if file exists). Falls back across multiple sources.
- `CellOntologyMapper(...)` defers heavy work: parsing OBO + computing embeddings happen on first `map_*` call. Pass `embeddings_path` to skip re-encode on subsequent runs.
- `model_name` defaults to `'all-mpnet-base-v2'` (~110M params); `'all-MiniLM-L6-v2'` is faster but ~5 % less accurate.
- `local_model_dir` is the HuggingFace cache path; set explicitly for offline / air-gapped use.

**Mapping algorithm:**
- Encode each input cell name with the sentence transformer.
- Cosine-similarity match against the encoded CL terms (and Cell Taxonomy entries when `use_taxonomy=True`).
- Pick the top-1 if `use_llm_selection=False`, OR pass top-`llm_candidates_count` to the LLM and let it pick when `use_llm_selection=True`.
- Below `threshold` → flag as `'No clear match'`.

**LLM expansion (`map_*_with_expansion`):**
- Called *before* embedding. The LLM rewrites short / abbreviation labels into descriptive text using `tissue_context` / `species` / `study_context` as prompt context.
- Cache: `setup_llm_expansion(cache_file=...)` JSON-persists prior expansions; same input string is never sent to the LLM twice.
- `api_type='custom_openai'` with `base_url=...` supports any OpenAI-compatible endpoint (Azure, Ollama, OhMyGPT, vLLM).

**Taxonomy mode:**
- `load_cell_taxonomy_resource(...)` reads the Jin 2023 Cell Taxonomy TSV; `species_filter` keeps only the relevant taxa to bound memory.
- `map_*_with_taxonomy(...)` combines CL match + Cell Taxonomy match; the `enhanced_*` obs columns expose both.

**Output obs columns** (after `map_adata`):
- `<new_col_name>` — best-match name (default `'cell_ontology'`).
- `<new_col_name>_cl_id` — CL ID (`CL:0000084`).
- `<new_col_name>_score` — cosine similarity.
- `<new_col_name>_ontology_id` — alternate ID (PopV-style if `load_ontology_mappings` was called).
- For `map_adata_with_taxonomy`: extra columns `_taxonomy_match` (CT term name) and `_ct_id` (CT ID).

## Notebook ↔ skill alignment

| Notebook section | Skill section |
|---|---|
| `download_cl(output_dir='new_ontology', filename='cl.json')` | Quick Workflow §1 |
| `CellOntologyMapper(cl_obo_file=..., embeddings_path=..., local_model_dir=...)` | Quick Workflow §2 |
| `mapper.map_adata(adata, cell_name_col='cell_label')` | Quick Workflow §3 (Mode 1) |
| `mapper.setup_llm_expansion(api_type='openai' / 'custom_openai', tissue_context='gut', species='mouse', ...)` | Quick Workflow §6 (Mode 2) |
| `mapper.map_adata_with_expansion(threshold=0.5, expand_abbreviations=True)` | Quick Workflow §7 |
| `mapper.load_cell_taxonomy_resource('Cell_Taxonomy_resource.txt', species_filter=...)` | Quick Workflow §10 |
| `mapper.map_adata_with_taxonomy(species='Mus musculus', tissue_context='Gut', threshold=0.3, ...)` | Quick Workflow §11 |
| `ov.pl.embedding(color=['cell_label', 'cell_ontology', 'enhanced_*'])` | Quick Workflow §12 (visualisation) |

## Docstring supplementation log

| Symbol | Prior state | Action |
|---|---|---|
| `CellOntologyMapper` (class) | 1-line ("🧬 Cell ontology mapping class using NLP") | filled — full Numpy-style docstring covering: NLP-encoder approach, three modes (plain / expansion / taxonomy), lifecycle (lazy embedding compute + persistence), output obs columns, typical workflow example. |

Method-level docstrings: most are 13–30 lines (good); a handful are 1-line (`clear_abbreviation_cache`, `load_embeddings`, `save_embeddings`, `save_mapping_results`, `get_statistics`, `get_ontology_statistics`, `print_mapping_summary*`, `show_expansion_summary`). Left as-is — these are simple persistence / inspection helpers whose names are self-explanatory and whose semantics are documented at the class level.

## Reviewer-Run Empirical Checks

- All cited functions importable: `from omicverse.single import CellOntologyMapper, download_cl` ✓
- The 28 public methods listed match `inspect.getmembers(CellOntologyMapper, predicate=inspect.isfunction)`.
- `setup_llm_expansion` `api_type` values verified: `'openai'` and `'custom_openai'` are the two paths in source.
- Output obs columns verified by reading the source: `<new_col_name>`, `<new_col_name>_cl_id`, `<new_col_name>_score`, `<new_col_name>_ontology_id` for plain; additional `_taxonomy_match` / `_ct_id` for taxonomy mode.
- No live smoke run executed; LLM expansion would require an OpenAI key.
