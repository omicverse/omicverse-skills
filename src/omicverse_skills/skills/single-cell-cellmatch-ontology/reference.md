# Cell Ontology mapping (CellMatch) — quick commands

## Plain mode

```python
import omicverse as ov
import scanpy as sc
ov.plot_set()

# 1) Download CL (idempotent, ~10 MB)
ov.single.download_cl(output_dir='new_ontology', filename='cl.json')

# 2) Construct mapper (first call computes ~3000-term embeddings, then cached)
mapper = ov.single.CellOntologyMapper(
    cl_obo_file='new_ontology/cl.json',
    embeddings_path='new_ontology/ontology_embeddings.pkl',
    local_model_dir='./my_models',
)

# 3) Map adata
mapping_results = mapper.map_adata(adata, cell_name_col='cell_label')
mapper.print_mapping_summary(mapping_results, top_n=15)

# 4) Inspect side-by-side
ov.pl.embedding(
    adata, basis='X_umap',
    color=['cell_label', 'cell_ontology', 'cell_ontology_cl_id'],
    wspace=0.55, ncols=2,
)
```

## LLM-expansion mode (for abbreviation-heavy labels)

```python
# OpenAI direct
mapper.setup_llm_expansion(
    api_type='openai',
    model='gpt-4o-2024-11-20',
    api_key='sk-...',
    tissue_context='gut',
    species='mouse',
    study_context='Epithelial cells from small intestine and organoids; some Salmonella or H. polygyrus infected.',
)

# OR custom OpenAI-compatible endpoint
mapper.setup_llm_expansion(
    api_type='custom_openai',
    api_key='sk-...',
    model='gpt-4.1-2025-04-14',
    base_url='https://api.ohmygpt.com/v1',
)

mapping_results = mapper.map_adata_with_expansion(
    adata=adata,
    cell_name_col='cell_label',
    threshold=0.5,
    expand_abbreviations=True,
)
mapper.print_mapping_summary(mapping_results, top_n=15)
mapper.show_expansion_summary(mapping_results)
```

## Taxonomy-enhanced mode (species + tissue)

```python
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

## Ad-hoc lookups

```python
# Find ontology terms most similar to a single name
top10 = mapper.find_similar_cells('TIL-1', top_k=10)
print(top10)

# Search by keyword
results = mapper.search_ontology_cells('memory T cell',
                                       case_sensitive=False, max_results=20)

# Browse by category
mapper.browse_ontology_by_category(categories=['T cell', 'B cell'],
                                   max_per_category=10)

# Marker-gene-based search (Cell Taxonomy)
ct_hits = mapper.search_by_marker(['CD8A', 'CD8B', 'GZMK'],
                                  species='Homo sapiens', top_k=10)

# Diagnostics
mapper.check_ontology_status()
print(mapper.get_ontology_statistics())
```

## Map a list of names without an AnnData

```python
names = ['TIL-1', 'TA-Early', 'IEL', 'plasmacytoid DC']
results = mapper.map_cells(names, threshold=0.5)
print(results)
```

## Persistence

```python
mapper.save_embeddings('cached_embeddings.pkl')
mapper.save_mapping_results(mapping_results, 'mapping_out.json')
mapper.clear_abbreviation_cache()      # next LLM run rebuilds
```
