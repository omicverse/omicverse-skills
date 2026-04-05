# Workflow Selection

## Reference build

- Use `workflow='standard'` for the notebook's default path.
- Use `nucleus=True` or `workflow='nucleus'` when the reference bundle is built for nucleus or velocity workflows.
- `workflow='lamanno'`, `workflow='nac'`, `workflow='kite'`, and `workflow='custom'` are source-present branches for specialized bundles.

## Counting

- Match `technology` to the chemistry of the FASTQs.
- Use `10XV3` for the notebook-like 10x V3 path.
- `BULK` bypasses barcode filtering and should not receive `filter_barcodes` or `filter_threshold`.
- Prefer `h5ad=True` for a quick downstream handoff when you do not need only text matrices.

## Flags to treat carefully

- `tcc=True` switches to transcript-compatibility counts.
- `mm=True` is disabled when `tcc=True`.
- `filter_barcodes=True` only makes sense for barcode-based single-cell libraries.
