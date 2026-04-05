from __future__ import annotations

import json
from pathlib import Path

import omicverse as ov
import scanpy as sc


REPO_ROOT = Path(__file__).resolve().parents[4]
DATA_PATH = REPO_ROOT / "data" / "hematopoiesis.h5ad"
WORK_ROOT = REPO_ROOT / "temp" / "cpdb_acceptance"
DB_PATH = WORK_ROOT / "cellphonedb.zip"
OUT_DIR = WORK_ROOT / "out"
TMP_DIR = WORK_ROOT / "tmp"


def main() -> int:
    WORK_ROOT.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(DATA_PATH)
    keep = adata.obs["cell_type"].value_counts().head(3).index.tolist()
    parts = [adata[adata.obs["cell_type"] == ct][:80] for ct in keep]
    adata = sc.concat(parts, join="inner")
    adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")

    cpdb_results, adata_cpdb = ov.single.run_cellphonedb_v5(
        adata,
        cpdb_file_path=str(DB_PATH),
        celltype_key="cell_type",
        min_cell_fraction=0.0,
        min_genes=50,
        min_cells=3,
        iterations=100,
        threshold=0.1,
        pvalue=0.1,
        threads=1,
        output_dir=str(OUT_DIR),
        temp_dir=str(TMP_DIR),
        cleanup_temp=False,
        debug=False,
    )

    viz = ov.pl.CellChatViz(adata_cpdb)
    count_matrix, weight_matrix = viz.compute_aggregated_network(
        pvalue_threshold=0.1,
        use_means=True,
    )
    pathway_comm = viz.compute_pathway_communication(
        method="mean",
        min_lr_pairs=1,
        min_expression=0.0,
    )
    _, _ = viz.get_significant_pathways_v2(
        pathway_comm,
        strength_threshold=0.0,
        pvalue_threshold=0.1,
        min_significant_pairs=1,
    )
    centrality = viz.netAnalysis_computeCentrality(
        signaling=None,
        pvalue_threshold=0.1,
        use_weight=True,
    )

    payload = {
        "count_shape": list(count_matrix.shape),
        "cpdb_obs": int(adata_cpdb.n_obs),
        "cpdb_vars": int(adata_cpdb.n_vars),
        "has_means_layer": "means" in adata_cpdb.layers,
        "has_pvalues_layer": "pvalues" in adata_cpdb.layers,
        "mode": "cellphonedb_v5",
        "n_pathways": len(pathway_comm),
        "sender_receiver_cols": all(
            key in adata_cpdb.obs.columns for key in ["sender", "receiver"]
        ),
        "centrality_keys": sorted(centrality.keys()),
        "weight_shape": list(weight_matrix.shape),
        "result_keys_ok": "means" in cpdb_results and "pvalues" in cpdb_results,
    }
    print(json.dumps(payload, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
