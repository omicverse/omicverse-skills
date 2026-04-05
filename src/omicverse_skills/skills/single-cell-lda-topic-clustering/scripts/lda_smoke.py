from __future__ import annotations

import json

import anndata as ad
import numpy as np
import pandas as pd

import omicverse as ov


def main() -> None:
    rng = np.random.default_rng(2)
    x = rng.poisson(4, size=(40, 25)).astype(float)
    x[:20, :5] += 10
    x[20:, 5:10] += 10
    adata = ad.AnnData(
        X=x,
        obs=pd.DataFrame(index=[f"c{i}" for i in range(40)]),
        var=pd.DataFrame(index=[f"g{i}" for i in range(25)]),
    )
    adata.var["highly_variable_features"] = True
    adata.layers["counts"] = adata.X.copy()
    ov.pp.normalize_total(adata)
    ov.pp.log1p(adata)
    ov.pp.scale(adata)
    ov.pp.pca(adata, n_pcs=8, layer="scaled")

    summary = {
        "callable_present": hasattr(ov.utils, "LDA_topic"),
        "counts_layer": "counts" in adata.layers,
        "has_pca": "X_pca" in adata.obsm,
    }
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
