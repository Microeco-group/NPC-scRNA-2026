# Script: 03_fusobacterium_distance.py
# Purpose: Distance from selected cell clusters to bacterial coordinates
# Workflow: spatial
# Required inputs: cell_anno_insitutype.h5ad; cosmx_tx_file.csv; optional insitutype_anno.csv
# Main outputs: output/tex_distance/disfuso_tex.h5ad; output/mreg_distance/disfuso_mreg.h5ad
# Prerequisites: Run from this workflow directory; configure input and output paths with command-line arguments.

import argparse
from pathlib import Path

import anndata as ad
import pandas as pd
from scipy.spatial import KDTree
import scanpy as sc


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute nearest-bacteria distances for selected cell clusters in a spatial AnnData object."
        )
    )
    parser.add_argument(
        "--adata",
        type=Path,
        default=Path("cell_anno_insitutype.h5ad"),
        help="Path to the annotated spatial AnnData object.",
    )
    parser.add_argument(
        "--annotation-csv",
        type=Path,
        default=None,
        help=(
            "Optional CSV containing cell_barcode and cluster_assignment columns to merge into adata.obs."
        ),
    )
    parser.add_argument(
        "--bacteria-csv",
        type=Path,
        default=Path("cosmx_tx_file.csv"),
        help="Path to the bacterial coordinate table.",
    )
    parser.add_argument(
        "--target",
        default="Fusobacterium genus",
        help="Target value in the bacteria table column 'target'.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("output"),
        help="Directory where output AnnData files will be written.",
    )
    parser.add_argument(
        "--tex-cluster",
        default="c07_CD8_Tex_PDCD1",
        help="Cluster label used for the Tex subset.",
    )
    parser.add_argument(
        "--mreg-cluster",
        default="c41_mregDC_LAMP3",
        help="Cluster label used for the mregDC subset.",
    )
    return parser.parse_args()


def add_cluster_annotations(adata_obj: ad.AnnData, annotation_csv: Path | None) -> ad.AnnData:
    if annotation_csv is None:
        if "cluster_assignment" not in adata_obj.obs.columns:
            raise ValueError(
                "cluster_assignment is missing from adata.obs and --annotation-csv was not provided."
            )
        return adata_obj

    anno = pd.read_csv(annotation_csv)
    required = {"cell_barcode", "cluster_assignment"}
    missing = required.difference(anno.columns)
    if missing:
        raise ValueError(f"Annotation CSV is missing required columns: {sorted(missing)}")

    anno_indexed = anno.set_index("cell_barcode")
    adata_obj.obs["cluster_assignment"] = anno_indexed.loc[adata_obj.obs_names, "cluster_assignment"].values
    return adata_obj


def compute_and_write_distances(
    adata_obj: ad.AnnData,
    bacteria_table: pd.DataFrame,
    cluster_label: str,
    output_path: Path,
) -> None:
    subset = adata_obj[adata_obj.obs["cluster_assignment"] == cluster_label].copy()
    if subset.n_obs == 0:
        raise ValueError(f"No cells matched cluster label: {cluster_label}")

    cell_coords = subset.obs[["CenterX_global_px", "CenterY_global_px"]].values
    bacteria_coords = bacteria_table[["x_global_px", "y_global_px"]].values
    if len(bacteria_coords) == 0:
        raise ValueError("No bacterial coordinates remain after filtering by target.")

    distances, _ = KDTree(bacteria_coords).query(cell_coords, k=1)
    subset.obs["dis2fuso"] = distances
    output_path.parent.mkdir(parents=True, exist_ok=True)
    subset.write(output_path)


def main() -> None:
    args = parse_args()
    adata_obj = sc.read_h5ad(args.adata)
    adata_obj = add_cluster_annotations(adata_obj, args.annotation_csv)

    bacteria = pd.read_csv(args.bacteria_csv, index_col=0)
    required_columns = {"target", "x_global_px", "y_global_px"}
    missing = required_columns.difference(bacteria.columns)
    if missing:
        raise ValueError(f"Bacteria table is missing required columns: {sorted(missing)}")

    bacteria = bacteria[bacteria["target"] == args.target].copy()
    compute_and_write_distances(
        adata_obj,
        bacteria,
        args.tex_cluster,
        args.output_root / "tex_distance" / "disfuso_tex.h5ad",
    )
    compute_and_write_distances(
        adata_obj,
        bacteria,
        args.mreg_cluster,
        args.output_root / "mreg_distance" / "disfuso_mreg.h5ad",
    )


if __name__ == "__main__":
    main()
