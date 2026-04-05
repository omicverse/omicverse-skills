# Source Grounding

## Inspected callables

- `ov.alignment.single.ref(...)`
- `ov.alignment.single.count(...)`
- `ov.alignment.ref(...)` and `ov.alignment.count(...)` as aliases
- `ov.alignment.single` namespace, built from `types.SimpleNamespace()`

## Source behavior worth preserving

- `ref` defaults to `workflow='standard'` and flips to `nucleus` when `nucleus=True` is passed with the default workflow.
- `ref` accepts several source-controlled workflow branches: `standard`, `nucleus`, `lamanno`, `nac`, `kite`, and `custom`.
- `count` defaults to `workflow='standard'` and adds `--filter bustools` only when `filter_barcodes=True` for non-BULK technologies.
- `count` raises if `technology='BULK'` is combined with `filter_barcodes=True` or `filter_threshold`.
- `count` uses `--mm` only when `tcc=False`.
- `count` returns discovered output file paths only if those files exist on disk.
- Both wrappers shell out to an external `kb` executable, falling back to `python -m kb` or `python -m kb_python` when available.

## Signature facts

- `ref(index_path, t2g_path, ..., workflow='standard', nucleus=False, ...)`
- `count(index_path, t2g_path, technology, fastq_paths, ..., h5ad=False, filter_barcodes=False, tcc=False, mm=False, nucleus=False, ...)`
