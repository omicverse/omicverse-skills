from __future__ import annotations

import json
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import omicverse as ov


def main() -> None:
    rng = np.random.default_rng(1)
    x = rng.poisson(lam=4, size=(120, 40)).astype(float)
    x[:60, :6] += 8
    x[60:, 6:12] += 8

    obs = pd.DataFrame(
        {
            "condition": ["Control"] * 60 + ["Salmonella"] * 60,
            "batch": ["b1"] * 20 + ["b2"] * 20 + ["b3"] * 20 + ["b4"] * 20 + ["b5"] * 20 + ["b6"] * 20,
            "cell_label": ["TA"] * 30 + ["Enterocyte"] * 30 + ["TA"] * 30 + ["Enterocyte"] * 30,
        },
        index=[f"cell{i}" for i in range(120)],
    )
    var = pd.DataFrame(index=[f"G{i}" for i in range(40)])
    adata = ad.AnnData(X=x, obs=obs, var=var)

    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.pca(adata, n_comps=10)
    adata.obsm["X_pca"] = adata.obsm["X_pca"].copy()

    dct = ov.single.DCT(
        adata,
        condition="condition",
        ctrl_group="Control",
        test_group="Salmonella",
        cell_type_key="cell_label",
        method="milopy",
        sample_key="batch",
        use_rep="X_pca",
    )

    summary = {
        "constructor_ok": True,
        "method": dct.method,
        "has_mdata": hasattr(dct, "mdata"),
        "has_milo_modality": hasattr(dct, "mdata") and "milo" in dct.mdata.mod,
        "conditions_used": sorted(dct.adata.obs["condition"].astype(str).unique().tolist()),
    }
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
