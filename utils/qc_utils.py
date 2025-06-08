import scanpy as sc
import anndata as an
import pandas as pd
import numpy as np

tresholds = {
    "total_counts": 1,
    "n_genes_by_counts": 1,
    "pct_counts_mt": 25,
    "pct_counts_ribo": 22,
    "pct_counts_hb": 2,
}

layer_palette = {
    "Empty spots": "#1f77b4",
    "WM": "#bcbd22",
    "L1": "#ff7f0e",
    "L6": "#e377c2",
    "L2": "#2ca02c",
    "L3": "#d62728",
    "L4": "#9467bd",
    "L5": "#8c564b",
    "L6a": "#ffc0cb",
    "L6b": "#fa9284",
}


def apply_upper_treshold(adata, column, treshold=None):
    global tresholds

    qc_name = "qc_" + column
    adata.obs[qc_name] = "True"
    treshold = tresholds[column] if treshold is None else treshold
    tresh_index = adata.obs[column] > treshold
    adata.obs.loc[tresh_index, qc_name] = "False"
    print(f"Threshold : {treshold}")
    print(f"Spots filtered: {tresh_index.sum()}")
    if tresh_index.sum():
        vc = adata.obs.groupby(qc_name)["label"].value_counts()
        print(vc[vc > 0]["False"])


def apply_percentile_treshold(adata, column, right_bound=True):
    global tresholds

    qc_name = "qc_" + column
    adata.obs[qc_name] = "True"
    tresh_num = 0

    for layer in adata.obs.label.unique():
        left_treshold = tresholds[column] / 100
        right_treshold = 1 - left_treshold if right_bound else 1
        left_bound = adata.obs[adata.obs.label == layer][column].quantile(left_treshold)
        right_bound = adata.obs[adata.obs.label == layer][column].quantile(
            right_treshold
        )

        tresh_index = (adata.obs.label == layer) & (
            (adata.obs[column] < left_bound) | (adata.obs[column] > right_bound)
        )
        adata.obs.loc[tresh_index, qc_name] = "False"
        tresh_num += tresh_index.sum()

    print(f"Spots filtered: {tresh_num}")
