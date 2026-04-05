from __future__ import annotations

import json
import os

import anndata as ad
import numpy as np
from scipy import sparse

import omicverse as ov


def build_smoke_adata() -> ad.AnnData:
    rng = np.random.default_rng(7)
    genes = [
        "Irf6",
        "Ets1",
        "E2f8",
        "Gata2",
        "Spi1",
        "Runx1",
        "Myc",
        "Fos",
        "Jun",
        "Stat1",
        "Klf1",
        "Tal1",
        "Lyl1",
        "Kit",
        "Gfi1b",
        "Hlf",
        "Mpl",
        "Elane",
        "Lmo2",
        "Cebpa",
    ]

    base = rng.poisson(lam=4, size=(24, len(genes))).astype(np.int64)
    base[:12, :6] += rng.poisson(lam=3, size=(12, 6))
    base[12:, 6:12] += rng.poisson(lam=3, size=(12, 6))
    counts = sparse.csr_matrix(base)

    adata = ad.AnnData(X=counts.copy())
    adata.var_names = genes
    adata.obs_names = [f"cell_{i}" for i in range(adata.n_obs)]
    adata.layers["counts"] = counts.copy()
    adata.layers["raw_count"] = counts.copy()
    adata.obs["cell_type"] = ["group_a"] * 12 + ["group_b"] * 12
    return adata


def main() -> int:
    db_glob = os.environ["SCENIC_DB_GLOB"]
    motif_path = os.environ["SCENIC_MOTIF_PATH"]

    adata = build_smoke_adata()
    scenic = ov.single.SCENIC(
        adata=adata,
        db_glob=db_glob,
        motif_path=motif_path,
        n_jobs=1,
    )
    edgelist = scenic.cal_grn(
        method="regdiffusion",
        layer="counts",
        n_steps=1,
        batch_size=8,
        device="cpu",
    )

    payload = {
        "db_glob_set": bool(db_glob),
        "motif_path_set": bool(motif_path),
        "cell_count": int(adata.n_obs),
        "gene_count": int(adata.n_vars),
        "edge_count": int(len(edgelist)),
        "columns_ok": list(edgelist.columns) == ["TF", "target", "importance"],
        "adjacency_synced": hasattr(scenic, "adjacencies") and len(scenic.adjacencies) == len(edgelist),
        "method": "regdiffusion",
        "layer": "counts",
    }
    print(json.dumps(payload, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
