from __future__ import annotations

import json
import os

import anndata as ad
import numpy as np
import pandas as pd

import omicverse as ov

os.environ.setdefault("JAX_PLATFORMS", "cpu")

from scib_metrics.benchmark import BatchCorrection, Benchmarker, BioConservation


def build_toy_adata() -> ad.AnnData:
    rng = np.random.default_rng(7)
    x = rng.poisson(4, size=(48, 36)).astype(float)
    x[:24, :6] += 5
    x[24:, 6:12] += 5
    x[::2, 12:18] += 3

    obs = pd.DataFrame(
        {
            "batch": pd.Categorical(["b0"] * 24 + ["b1"] * 24),
            "cell_type": pd.Categorical(["t0"] * 12 + ["t1"] * 12 + ["t0"] * 12 + ["t1"] * 12),
        },
        index=[f"c{i}" for i in range(48)],
    )
    var = pd.DataFrame(index=[f"g{i}" for i in range(36)])
    adata = ad.AnnData(X=x, obs=obs, var=var)
    adata.layers["counts"] = adata.X.copy()
    adata.var["highly_variable_features"] = True

    ov.pp.normalize_total(adata)
    ov.pp.log1p(adata)
    ov.pp.scale(adata)
    ov.pp.pca(adata, n_pcs=8, layer="scaled")
    adata.obsm["X_pca"] = adata.obsm["scaled|original|X_pca"].copy()
    return adata


def main() -> None:
    adata = build_toy_adata()

    ov.single.batch_correction(
        adata,
        batch_key="batch",
        methods="harmony",
        n_pcs=8,
        use_gpu=False,
        verbose=False,
    )
    ov.single.batch_correction(
        adata,
        batch_key="batch",
        methods="combat",
        n_pcs=8,
    )

    bm = Benchmarker(
        adata,
        batch_key="batch",
        label_key="cell_type",
        embedding_obsm_keys=["X_pca", "X_pca_harmony", "X_combat"],
        bio_conservation_metrics=BioConservation(
            isolated_labels=False,
            nmi_ari_cluster_labels_leiden=False,
            nmi_ari_cluster_labels_kmeans=False,
            silhouette_label=True,
            clisi_knn=False,
        ),
        batch_correction_metrics=BatchCorrection(
            silhouette_batch=True,
            ilisi_knn=False,
            kbet_per_label=False,
            graph_connectivity=False,
            pcr_comparison=True,
        ),
        pre_integrated_embedding_obsm_key="X_pca",
        n_jobs=1,
        progress_bar=False,
    )
    bm.benchmark()
    results = bm.get_results(min_max_scale=False)

    summary = {
        "x_pca": "X_pca" in adata.obsm,
        "x_pca_harmony": "X_pca_harmony" in adata.obsm,
        "x_harmony": "X_harmony" in adata.obsm,
        "x_combat": "X_combat" in adata.obsm,
        "benchmark_nonempty": results.shape[0] > 0,
        "benchmark_has_harmony": "X_pca_harmony" in results.index,
        "benchmark_has_combat": "X_combat" in results.index,
    }
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
