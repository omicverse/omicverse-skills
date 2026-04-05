from __future__ import annotations

import json
import numpy as np
import pandas as pd
import anndata as ad
import omicverse as ov


def main() -> None:
    rng = np.random.default_rng(0)
    x = rng.poisson(lam=4, size=(80, 30)).astype(float)
    x[:40, :5] += 12
    x[40:, 5:10] += 12

    obs = pd.DataFrame(
        {
            "condition": ["Control"] * 40 + ["Salmonella"] * 40,
            "cell_label": ["TA"] * 20 + ["Enterocyte"] * 20 + ["TA"] * 20 + ["Enterocyte"] * 20,
        },
        index=[f"cell{i}" for i in range(80)],
    )
    var = pd.DataFrame(index=[f"G{i}" for i in range(30)])
    adata = ad.AnnData(X=x, obs=obs, var=var)

    deg = ov.single.DEG(
        adata,
        condition="condition",
        ctrl_group="Control",
        test_group="Salmonella",
        method="wilcoxon",
        use_raw=False,
    )
    deg.run(celltype_key="cell_label", celltype_group=["TA"], max_cells=1000)
    res = deg.get_results()

    summary = {
        "rows": int(res.shape[0]),
        "has_required_columns": all(
            key in res.columns
            for key in ["log2FC", "pvalue", "padj", "baseMean", "pct_ctrl", "pct_test", "pct_diff"]
        ),
        "celltypes_used": sorted(deg.adata_test.obs["cell_label"].astype(str).unique().tolist()),
        "conditions_used": sorted(deg.adata_test.obs["condition"].astype(str).unique().tolist()),
    }
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
