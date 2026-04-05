from __future__ import annotations

import argparse
import json

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

import omicverse as ov


def make_fixture() -> ad.AnnData:
    rng = np.random.default_rng(0)
    clusters = np.repeat(["Ductal", "Alpha", "Beta", "Delta"], 20)
    x = rng.poisson(lam=3, size=(len(clusters), 60)).astype(float)
    x[clusters == "Ductal", :8] += 20
    x[clusters == "Alpha", 8:16] += 20
    x[clusters == "Beta", 16:24] += 20
    x[clusters == "Delta", 24:32] += 20

    adata = ad.AnnData(
        X=x,
        obs=pd.DataFrame(
            {"clusters": pd.Categorical(clusters)},
            index=[f"cell{i}" for i in range(len(clusters))],
        ),
        var=pd.DataFrame(index=[f"G{i}" for i in range(60)]),
    )
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=12)
    adata.obsm["X_umap"] = adata.obsm["X_pca"][:, :2].copy()
    sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=10)
    return adata


def run_diffusion_map() -> dict[str, object]:
    adata = make_fixture()
    traj = ov.single.TrajInfer(
        adata,
        basis="X_umap",
        use_rep="X_pca",
        n_comps=12,
        n_neighbors=10,
        groupby="clusters",
    )
    traj.set_origin_cells("Ductal")
    traj.inference(method="diffusion_map")
    ov.utils.cal_paga(
        adata,
        groups="clusters",
        vkey="paga",
        use_time_prior="dpt_pseudotime",
    )
    return {
        "branch": "diffusion_map",
        "has_dpt_pseudotime": "dpt_pseudotime" in adata.obs,
        "has_paga": "paga" in adata.uns,
        "paga_groups": adata.uns["paga"]["groups"],
    }


def run_palantir() -> dict[str, object]:
    import omicverse.external.palantir.utils as pal_utils

    orig_magic = pal_utils.run_magic_imputation

    def patched_run_magic_imputation(*args, **kwargs):
        kwargs.setdefault("n_jobs", 1)
        return orig_magic(*args, **kwargs)

    pal_utils.run_magic_imputation = patched_run_magic_imputation

    adata = make_fixture()
    traj = ov.single.TrajInfer(
        adata,
        basis="X_umap",
        use_rep="X_pca",
        n_comps=12,
        n_neighbors=10,
        groupby="clusters",
    )
    traj.set_origin_cells("Ductal")
    traj.set_terminal_cells(["Alpha", "Beta", "Delta"])
    traj.inference(method="palantir", num_waypoints=40, knn=10, n_jobs=1)
    traj.palantir_cal_branch(eps=0)
    ov.utils.cal_paga(
        adata,
        groups="clusters",
        vkey="paga",
        use_time_prior="palantir_pseudotime",
    )
    return {
        "branch": "palantir",
        "patched_magic_single_process": True,
        "has_palantir_pseudotime": "palantir_pseudotime" in adata.obs,
        "has_palantir_entropy": "palantir_entropy" in adata.obs,
        "has_fate_probabilities": "palantir_fate_probabilities" in adata.obsm,
        "has_branch_masks": "branch_masks" in adata.obsm,
        "has_paga": "paga" in adata.uns,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--branch", choices=("diffusion_map", "palantir"), required=True)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.branch == "diffusion_map":
        payload = run_diffusion_map()
    else:
        payload = run_palantir()
    print(json.dumps(payload, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
