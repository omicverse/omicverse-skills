from __future__ import annotations

import json
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

import omicverse as ov


def main() -> None:
    rng = np.random.default_rng(3)
    x = rng.poisson(4, size=(40, 30)).astype(float)
    x[:20, :5] += 8
    x[20:, 5:10] += 8
    adata = ad.AnnData(
        X=x,
        obs=pd.DataFrame(index=[f"c{i}" for i in range(40)]),
        var=pd.DataFrame(index=[f"g{i}" for i in range(30)]),
    )
    adata.layers["counts"] = adata.X.copy()
    ov.pp.normalize_total(adata)
    ov.pp.log1p(adata)
    ov.pp.scale(adata)
    ov.pp.pca(adata, n_pcs=8, layer="scaled")

    outdir = Path(tempfile.mkdtemp(prefix="cnmf_smoke_"))
    cnmf = ov.single.cNMF(
        adata,
        components=np.array([2]),
        n_iter=1,
        seed=14,
        num_highvar_genes=20,
        output_dir=str(outdir),
        name="toy",
        use_gpu=False,
    )
    summary = {
        "constructed": True,
        "use_gpu": cnmf.use_gpu,
        "name": cnmf.name,
    }
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
