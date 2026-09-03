#!/usr/bin/env python3
"""Memory-conscious exploratory analysis of the GSE72056 melanoma scRNA-seq matrix."""

from __future__ import annotations

import argparse
import csv
import gzip
import importlib.metadata
import json
import os
import tempfile
import urllib.request
import warnings
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "alex_w0731_mpl_cache"))
os.environ.setdefault("LOKY_MAX_CPU_COUNT", "4")

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests


DATA_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/GSE72056/suppl/"
    "GSE72056_melanoma_single_cell_revised_v2.txt.gz"
)
META_TUMOR = "tumor"
META_MALIGNANT = "malignant(1=no,2=yes,0=unresolved)"
META_CELL_TYPE = "non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)"

MARKER_GENES = {
    "T/NK": ["CD3D", "CD3E", "TRAC", "NKG7", "GNLY"],
    "B": ["CD79A", "MS4A1", "CD37", "CD74"],
    "Myeloid": ["LST1", "TYROBP", "FCER1G", "CTSS"],
    "Endothelial": ["PECAM1", "VWF", "KDR", "EMCN"],
    "CAF": ["COL1A1", "COL1A2", "DCN", "COL3A1"],
    "Melanoma": ["MLANA", "PMEL", "MITF", "TYR"],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("data/GSE72056_melanoma_single_cell_revised_v2.txt.gz"),
        help="Path to the GEO matrix (downloaded automatically when absent).",
    )
    parser.add_argument("--hvg", type=int, default=1500, help="Number of variable genes to retain.")
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args()


def ensure_input(path: Path) -> None:
    if path.exists():
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading public GEO matrix to {path} ...")
    urllib.request.urlretrieve(DATA_URL, path)


def open_rows(path: Path):
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        yield "__HEADER__", header[1:]
        for row in reader:
            if len(row) < 2:
                continue
            yield row[0].strip('"'), row[1:]


def numeric_vector(values: list[str], n_cells: int) -> np.ndarray:
    vector = np.fromiter((float(value) if value else 0.0 for value in values), dtype=np.float32)
    if vector.size != n_cells:
        raise ValueError(f"Expected {n_cells} values, observed {vector.size}.")
    return vector


def first_pass(path: Path, hvg_count: int):
    iterator = open_rows(path)
    _, cell_ids = next(iterator)
    n_cells = len(cell_ids)
    metadata: dict[str, list[str]] = {}
    stats: list[tuple[str, float, float, float]] = []
    detected_genes = np.zeros(n_cells, dtype=np.int32)
    total_expression = np.zeros(n_cells, dtype=np.float64)

    for gene, values in iterator:
        if gene in {META_TUMOR, META_MALIGNANT, META_CELL_TYPE}:
            metadata[gene] = values
            continue

        vector = numeric_vector(values, n_cells)
        detected_fraction = float(np.mean(vector > 0))
        detected_genes += vector > 0
        total_expression += vector

        if detected_fraction >= 0.01 and not gene.startswith(("MT-", "RPL", "RPS")):
            stats.append((gene, float(vector.mean()), float(vector.var()), detected_fraction))

    missing_meta = {META_TUMOR, META_MALIGNANT, META_CELL_TYPE} - set(metadata)
    if missing_meta:
        raise ValueError(f"Required metadata rows are missing: {sorted(missing_meta)}")

    stats_frame = pd.DataFrame(stats, columns=["gene", "mean", "variance", "detected_fraction"])
    stats_frame = stats_frame.sort_values("variance", ascending=False)
    selected = set(stats_frame.head(hvg_count)["gene"])
    selected.update(gene for genes in MARKER_GENES.values() for gene in genes)
    return cell_ids, metadata, stats_frame, selected, detected_genes, total_expression


def second_pass(path: Path, selected: set[str], n_cells: int):
    genes: list[str] = []
    rows: list[np.ndarray] = []
    iterator = open_rows(path)
    next(iterator)
    for gene, values in iterator:
        if gene in selected:
            genes.append(gene)
            rows.append(numeric_vector(values, n_cells))
    if not rows:
        raise ValueError("No selected genes were found in the expression matrix.")
    return genes, np.vstack(rows).T


def build_annotations(metadata: dict[str, list[str]]) -> pd.DataFrame:
    malignant_map = {"0": "Unresolved", "1": "Non-malignant", "2": "Malignant"}
    lineage_map = {
        "0": "Unresolved",
        "1": "T cell",
        "2": "B cell",
        "3": "Macrophage",
        "4": "Endothelial",
        "5": "CAF",
        "6": "NK cell",
    }
    malignant = pd.Series(metadata[META_MALIGNANT]).map(malignant_map).fillna("Unresolved")
    lineage = pd.Series(metadata[META_CELL_TYPE]).map(lineage_map).fillna("Unresolved")
    annotation = np.where(malignant.eq("Malignant"), "Malignant", lineage)
    return pd.DataFrame(
        {
            "patient": [f"Patient {value}" for value in metadata[META_TUMOR]],
            "malignancy": malignant.to_numpy(),
            "annotation": annotation,
        }
    )


def save_umap(adata: ad.AnnData, color_by: str, path: Path, title: str) -> None:
    frame = pd.DataFrame(adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"], index=adata.obs_names)
    frame[color_by] = adata.obs[color_by].astype(str).to_numpy()
    categories = sorted(frame[color_by].unique())
    if categories and all(category.isdigit() for category in categories):
        categories = sorted(categories, key=int)
    palette = plt.get_cmap("tab20", max(len(categories), 1))

    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    for index, category in enumerate(categories):
        group = frame[frame[color_by] == category]
        ax.scatter(group.UMAP1, group.UMAP2, s=5, alpha=0.72, label=category, color=palette(index))
    ax.set_title(title, loc="left", weight="bold")
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left", markerscale=2)
    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def rank_cluster_markers(adata: ad.AnnData) -> pd.DataFrame:
    """Rank markers with Welch tests on the published log2-expression scale."""
    matrix = np.asarray(adata.raw.X)
    clusters = sorted(adata.obs["leiden"].astype(str).unique(), key=int)
    frames: list[pd.DataFrame] = []

    for cluster in clusters:
        inside = adata.obs["leiden"].astype(str).to_numpy() == cluster
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            statistic, p_value = ttest_ind(
                matrix[inside],
                matrix[~inside],
                axis=0,
                equal_var=False,
                nan_policy="omit",
            )
        p_value = np.where(np.isfinite(p_value), p_value, 1.0)
        p_value = np.clip(p_value, np.finfo(float).tiny, 1.0)
        frame = pd.DataFrame(
            {
                "cluster": cluster,
                "gene": adata.raw.var_names,
                "cells_in_cluster": int(inside.sum()),
                "mean_in_cluster": matrix[inside].mean(axis=0),
                "mean_outside_cluster": matrix[~inside].mean(axis=0),
                "t_statistic": statistic,
                "p_value": p_value,
                "fdr": multipletests(p_value, method="fdr_bh")[1],
            }
        )
        frame["mean_log2_difference"] = frame["mean_in_cluster"] - frame["mean_outside_cluster"]
        frames.append(frame.sort_values(["fdr", "p_value", "mean_log2_difference"], ascending=[True, True, False]).head(10))

    return pd.concat(frames, ignore_index=True)


def main() -> None:
    args = parse_args()
    np.random.seed(args.seed)
    sc.settings.seed = args.seed
    repo_root = Path(__file__).resolve().parents[1]
    figure_dir = repo_root / "results" / "figures"
    table_dir = repo_root / "results" / "tables"
    processed_dir = repo_root / "data" / "processed"
    for directory in (figure_dir, table_dir, processed_dir):
        directory.mkdir(parents=True, exist_ok=True)

    ensure_input(args.input)
    print("Pass 1/2: collecting QC metrics and ranking variable genes...")
    cell_ids, metadata, gene_stats, selected, detected, total = first_pass(args.input, args.hvg)
    print("Pass 2/2: constructing the selected expression matrix...")
    genes, matrix = second_pass(args.input, selected, len(cell_ids))

    obs = build_annotations(metadata)
    obs.index = pd.Index(cell_ids, name="cell_id")
    obs["n_genes_detected"] = detected
    obs["total_log_expression"] = total
    var = pd.DataFrame(index=pd.Index(genes, name="gene"))
    adata = ad.AnnData(X=matrix, obs=obs, var=var)
    adata.uns["log1p"] = {"base": 2}
    adata.raw = adata.copy()

    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=40, svd_solver="arpack", random_state=args.seed)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=args.seed)
    sc.tl.umap(adata, random_state=args.seed)
    sc.tl.leiden(
        adata,
        resolution=0.6,
        random_state=args.seed,
        key_added="leiden",
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )

    save_umap(
        adata,
        "annotation",
        figure_dir / "umap_author_annotation.png",
        "GSE72056 melanoma cells by published annotation",
    )
    save_umap(
        adata,
        "leiden",
        figure_dir / "umap_leiden_clusters.png",
        "Unsupervised Leiden clusters",
    )

    available_markers = {
        group: [gene for gene in markers if gene in adata.raw.var_names]
        for group, markers in MARKER_GENES.items()
    }
    available_markers = {group: genes for group, genes in available_markers.items() if genes}
    dotplot = sc.pl.dotplot(
        adata,
        available_markers,
        groupby="annotation",
        use_raw=True,
        show=False,
        return_fig=True,
        standard_scale="var",
        title="Canonical lineage-marker expression",
    )
    dotplot.savefig(figure_dir / "marker_dotplot.png", dpi=200)

    rank_cluster_markers(adata).to_csv(table_dir / "top_markers_by_cluster.csv", index=False)

    adata.obs["annotation"].value_counts().rename_axis("annotation").reset_index(name="cells").to_csv(
        table_dir / "cell_counts_by_annotation.csv", index=False
    )
    pd.crosstab(adata.obs["leiden"], adata.obs["annotation"]).to_csv(
        table_dir / "cluster_annotation_crosstab.csv"
    )
    gene_stats.head(args.hvg).to_csv(table_dir / "selected_variable_genes.csv", index=False)

    summary = {
        "accession": "GSE72056",
        "cells": int(adata.n_obs),
        "genes_in_working_matrix": int(adata.n_vars),
        "patients": int(adata.obs["patient"].nunique()),
        "leiden_clusters": int(adata.obs["leiden"].nunique()),
        "seed": args.seed,
        "input_scale": "Published processed log-expression matrix; no count normalization repeated",
    }
    (repo_root / "results" / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    versions = {
        package: importlib.metadata.version(package)
        for package in ("anndata", "igraph", "leidenalg", "matplotlib", "numpy", "pandas", "scanpy", "scipy", "statsmodels")
    }
    (repo_root / "results" / "software_versions.json").write_text(
        json.dumps(versions, indent=2) + "\n", encoding="utf-8"
    )
    adata.write_h5ad(processed_dir / "gse72056_hvg_processed.h5ad", compression="gzip")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
