from __future__ import annotations

import json

import anndata as ad
import numpy as np
import pandas as pd

import omicverse as ov


def main() -> None:
    rng = np.random.default_rng(0)
    x = rng.poisson(4, size=(80, 30)).astype(float)
    x[:40, :5] += 8
    x[40:, 5:10] += 8
    adata = ad.AnnData(
        X=x,
        obs=pd.DataFrame(index=[f"c{i}" for i in range(80)]),
        var=pd.DataFrame(index=[f"g{i}" for i in range(30)]),
    )

    ov.pp.normalize_total(adata)
    ov.pp.log1p(adata)
    ov.pp.scale(adata)
    ov.pp.pca(adata, n_pcs=10, layer="scaled")
    ov.pp.neighbors(adata, n_neighbors=10, n_pcs=10, use_rep="X_pca")

    ov.utils.cluster(adata, method="leiden", resolution=0.5)
    ov.utils.cluster(adata, method="GMM", use_rep="X_pca", n_components=3)

    summary = {
        "leiden": "leiden" in adata.obs.columns,
        "gmm_cluster": "gmm_cluster" in adata.obs.columns,
        "mclust": "mclust" in adata.obs.columns,
    }
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
