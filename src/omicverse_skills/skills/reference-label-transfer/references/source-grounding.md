# Source Grounding

This skill is grounded in the live OmicVerse source and inspected runtime signatures.

## Inspection Notes

- Source inspected from the local OmicVerse checkout.
- Interface inspection ran in `omictest`.
- The reference-transfer workflow was validated with a synthetic query/reference pair before writing the skill.

## Primary Source Files

- `omicverse/single/_annotation_ref.py`
- `omicverse/single/_batch.py`
- `omicverse/pp/_preprocess.py`
- `omicverse/pp/_qc.py`

## Inspected Signatures

### `AnnotationRef`

Inspected signature:

```text
(adata_query: anndata.AnnData, adata_ref: anndata.AnnData, celltype_key: str = 'celltype')
```

Source behavior:

- checks shared gene names between query and reference
- concatenates the two objects into `self.adata_new`
- warns when the objects appear to be log-normalized and should be recovered to counts

### `AnnotationRef.preprocess`

Inspected signature:

```text
(self, mode='shiftlog|pearson', n_HVGs=3000, batch_key='integrate_batch')
```

Notebook-covered branch:

- `shiftlog|pearson`

Source-level preprocessing mode options from `omicverse.pp.preprocess(...)`:

- `shiftlog|pearson`
- `shiftlog|seurat`
- `pearson|pearson`
- `pearson|seurat`

### `AnnotationRef.train`

Inspected signature:

```text
(self, method='harmony', **kwargs)
```

Source-grounded branch values:

- `harmony`
- `scVI`
- `scanorama`

Notebook-covered branches:

- `harmony`
- `scVI`

### `AnnotationRef.predict`

Inspected signature:

```text
(self, method='harmony', n_neighbors=15, pred_key=None, uncert_key=None)
```

Source-grounded branch values:

- `harmony`
- `scVI`
- `scanorama`

### `batch_correction`

Inspected signature:

```text
(adata, batch_key, use_rep='scaled|original|X_pca', methods='harmony', n_pcs=50, **kwargs)
```

Source-grounded `methods` branches:

- `harmony`
- `combat`
- `scanorama`
- `scVI`
- `CellANOVA`
- `Concord`

The notebook's reference-transfer workflow uses the `harmony` and `scVI` transfer branches only, but the underlying integration helper exposes more methods.

## Source-Derived Behavior Notes

### Harmony branch

- `train(method='harmony')` stores `X_pca_harmony_anno` on both query and reference objects
- `predict(method='harmony')` writes `harmony_prediction` and `harmony_uncertainty`
- the default `n_neighbors` is 15

### scVI branch

- `train(method='scVI')` stores `X_scVI_anno` on both query and reference objects
- `predict(method='scVI')` writes `scVI_prediction` and `scVI_uncertainty`

### Scanorama branch

- the live source exposes the branch even though the notebook does not make it its primary path
- the transfer keys become `scanorama_prediction` and `scanorama_uncertainty`

### Preprocessing handoff

- `AnnotationRef.preprocess(...)` forwards its `mode` to the shared preprocessing helper
- the notebook's `mode='shiftlog|pearson'` path preserves the raw-count + Pearson-residual style workflow

## Interface Inspection Evidence

The following live signatures were checked directly in `omictest`:

- `AnnotationRef`
- `AnnotationRef.preprocess`
- `AnnotationRef.train`
- `AnnotationRef.predict`

## Human Review Notes

- The notebook's reference-selection cells are important, but the reusable skill should start at the paired query/reference stage.
- The synthetic smoke path validated harmony transfer on a tiny paired dataset and confirmed that prediction columns and integrated embedding keys are produced.
