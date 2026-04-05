from __future__ import annotations

import json
import zipfile
from pathlib import Path

import omicverse as ov
import scanpy as sc


REPO_ROOT = Path(__file__).resolve().parents[4]
DATA_PATH = REPO_ROOT / "data" / "DentateGyrus" / "10X43_1.h5ad"
WORK_ROOT = REPO_ROOT / "temp" / "cytotrace2_acceptance"
ZIP_PATH = WORK_ROOT / "cymodels.zip"
MODEL_ROOT = WORK_ROOT / "models"
MODEL_DIR = MODEL_ROOT / "5_models_weights"
OUT_DIR = WORK_ROOT / "out"


def ensure_models() -> None:
    WORK_ROOT.mkdir(parents=True, exist_ok=True)
    if MODEL_DIR.exists():
        return
    if not ZIP_PATH.exists():
        ov.datasets.download_data(
            "https://stacks.stanford.edu/file/cv694yk7414/cymodels.zip",
            dir=str(WORK_ROOT),
        )
    with zipfile.ZipFile(ZIP_PATH, "r") as handle:
        handle.extractall(MODEL_ROOT)


def main() -> int:
    ensure_models()

    adata = sc.read_h5ad(DATA_PATH)
    adata = adata[:120, :].copy()
    adata = ov.pp.preprocess(adata, mode="shiftlog|pearson", n_HVGs=500)
    results = ov.single.cytotrace2(
        adata,
        use_model_dir=str(MODEL_DIR),
        species="mouse",
        batch_size=120,
        smooth_batch_size=60,
        disable_parallelization=True,
        max_pcs=50,
        output_dir=str(OUT_DIR),
    )

    payload = {
        "rows": int(results.shape[0]),
        "columns_ok": list(results.columns)
        == [
            "preKNN_CytoTRACE2_Score",
            "preKNN_CytoTRACE2_Potency",
            "CytoTRACE2_Score",
            "CytoTRACE2_Potency",
            "CytoTRACE2_Relative",
        ],
        "obs_keys_ok": all(
            key in adata.obs.columns
            for key in [
                "CytoTRACE2_Score",
                "CytoTRACE2_Potency",
                "CytoTRACE2_Relative",
                "preKNN_CytoTRACE2_Score",
                "preKNN_CytoTRACE2_Potency",
            ]
        ),
        "result_file": (OUT_DIR / "cytotrace2_results.txt").exists(),
        "mode": "shiftlog|pearson",
        "species": "mouse",
    }
    print(json.dumps(payload, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
