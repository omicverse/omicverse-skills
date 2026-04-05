from __future__ import annotations

import json
import multiprocessing as mp

mp.set_start_method("fork", force=True)

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import anndata as ad
import numpy as np

import omicverse as ov


def main() -> int:
    rng = np.random.default_rng(7)
    X = rng.poisson(3, size=(80, 120)).astype(float)
    adata = ad.AnnData(X=X)
    adata.var_names = [f"g{i}" for i in range(adata.n_vars)]
    adata.obs_names = [f"c{i}" for i in range(adata.n_obs)]
    adata.layers["spliced"] = rng.poisson(3, size=X.shape).astype(float)
    adata.layers["unspliced"] = rng.poisson(1, size=X.shape).astype(float)

    ov.settings.cpu_init()
    velo = ov.single.Velo(adata)
    velo.filter_genes(min_shared_counts=3)
    velo.preprocess(recipe="monocle", n_neighbors=15, n_pcs=20)
    velo.moments(backend="scvelo", n_neighbors=15, n_pcs=20)
    velo.dynamics(backend="scvelo", fit_scaling=False, n_jobs=1)
    velo.cal_velocity(method="scvelo")
    velo.velocity_graph(vkey="velocity_S", xkey="Ms", n_jobs=1)
    ov.pp.umap(adata)
    velo.velocity_embedding(basis="umap", vkey="velocity_S")

    # Use a broad subset so the GraphVelo PCA has enough features on synthetic data.
    gene_subset = velo.adata.var_names.tolist()
    velo.graphvelo(
        xkey="Ms",
        vkey="velocity_S",
        n_jobs=1,
        basis_keys=["X_umap", "X_pca"],
        gene_subset=gene_subset,
    )
    velo.velocity_graph(vkey="velocity_gv", xkey="Ms", n_jobs=1)
    velo.velocity_embedding(basis="umap", vkey="velocity_gv")

    fig, ax = plt.subplots(figsize=(4, 4))
    ov.pl.embedding(
        adata,
        basis="X_umap",
        color="leiden" if "leiden" in adata.obs else None,
        ax=ax,
        show=False,
        size=45,
        alpha=0.3,
    )
    try:
        ov.pl.add_streamplot(adata, basis="X_umap", velocity_key="velocity_gv_umap", ax=ax)
        stream_ok = True
    except Exception as exc:  # pragma: no cover - smoke only
        stream_ok = repr(exc)
    plt.close(fig)

    summary = {
        "x_pca": "X_pca" in adata.obsm,
        "x_umap": "X_umap" in adata.obsm,
        "velocity_s": "velocity_S" in adata.layers,
        "velocity_s_graph": "velocity_S_graph" in adata.uns,
        "velocity_s_umap": "velocity_S_umap" in adata.obsm,
        "velocity_gv": "velocity_gv" in adata.layers,
        "velocity_gv_genes": "velocity_gv_genes" in adata.var,
        "velocity_gv_graph": "velocity_gv_graph" in adata.uns,
        "velocity_gv_umap": "velocity_gv_umap" in adata.obsm,
        "gv_x_umap": "gv_X_umap" in adata.obsm,
        "gv_x_pca": "gv_X_pca" in adata.obsm,
        "streamplot_ok": stream_ok is True,
    }
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
