# Compatibility Notes

## Notebook Drift vs Reusable Contract

- The notebook demonstrates one dataset-specific subset and one cell type, but the skill generalizes this to any compatible condition labels and cell-type labels.
- The notebook presents memento as a second DEG method. In current code, that branch has a different data contract because it may need counts, a `counts` layer, or count recovery.

## Runtime Considerations

- Large subsets can trigger downsampling through `max_cells`; that is current wrapper behavior, not a notebook-only heuristic.
- The result columns for the scanpy-backed branches are added by OmicVerse after `rank_genes_groups`, so downstream code should expect the OmicVerse-formatted table rather than raw scanpy output.
- Keep acceptance and smoke execution shell-neutral. The generated skill only assumes a POSIX-like shell can invoke `${PYTHON} script.py`; it does not require `zsh` features.
