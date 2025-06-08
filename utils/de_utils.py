import itertools
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
from matplotlib import ticker
import matplotlib.pyplot as plt
from matplotlib.patches import bbox_artist
from statsmodels.stats.multitest import multipletests
from mpl_toolkits.axes_grid1 import make_axes_locatable

layers = ["L1", "L2", "L3", "L4", "L5", "L6", "WM"]


def heatmap_annot(data, annotation, fig, ax, title="", **kwargs):
    """
    Create a heatmap with annotations on top
    """
    sns.heatmap(data, ax=ax, **kwargs)
    ax.tick_params(axis="both", which="major", labelsize=10)
    ax.grid(False)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("top", size=0.15, pad=0.05)
    cmap = plt.get_cmap("Set3")

    annot_size = annotation.value_counts().loc[layers].values
    annot_pos = annot_size.cumsum()
    bounds = [0] + list(annot_pos)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    fig.colorbar(
        mpl.cm.ScalarMappable(cmap=cmap, norm=norm),
        cax=cax,
        ticks=bounds,
        orientation="horizontal",
        spacing="proportional",
    )
    cax.xaxis.set_major_locator(ticker.FixedLocator(annot_pos - annot_size / 2))
    cax.xaxis.set_major_formatter(ticker.FixedFormatter(layers))
    cax.xaxis.tick_top()

    cax.set_title(title)


def summary_de(de_dict):
    """
    Summarize the number of upregulated, downregulated, and not significant genes
    in each layer from the DE results dictionary.
    """
    summary_dict = dict()
    for layer, df in de_dict.items():
        up_reg = df[(df.p_val_adj < 0.05) & (df.logFC > 0)].shape[0]
        down_reg = df[(df.p_val_adj < 0.05) & (df.logFC < 0)].shape[0]
        not_sign = df.shape[0] - up_reg - down_reg
        summary_dict[layer] = [up_reg, not_sign, down_reg]
    summary_df = pd.DataFrame(summary_dict, index=["up_reg", "not_sign", "down_reg"])
    return summary_df


def count_unique(marker_dict):
    """
    Count unique genes in each layer and genes shared between layers.
    """
    countgene_list = []
    unique_genes = dict()
    for layer in layers:
        layer_list = []
        for oth_lay in layers:
            if layer == oth_lay:
                marker_genes_other = list(
                    itertools.chain(
                        *[
                            marker_dict[oth_lay]
                            for oth_lay in layers
                            if oth_lay != layer
                        ]
                    )
                )
                unique_genes[layer] = set(marker_dict[layer]).difference(
                    marker_genes_other
                )
                layer_list.append(len(unique_genes[layer]))
            else:
                layer_list.append(
                    len(set(marker_dict[layer]).intersection(marker_dict[oth_lay]))
                )
        countgene_list.append(layer_list)

    countgene_df = pd.DataFrame(countgene_list, index=layers, columns=layers)
    return countgene_df, unique_genes


def multipletest_correction(dict_with_df):
    """
    Apply multiple testing correction to the P-values in each DataFrame
    """
    for layer, df in dict_with_df.items():
        mult_test = multipletests(df["PValue"], method="fdr_bh")
        df["p_val_adj"] = mult_test[1]


def calculate_layerpct(adata, layer):
    """
    Calculate percent of spots where gene's expressed for the certain layer
    """
    adata_subs = adata[adata.obs.label == layer]
    n_obs = adata_subs.n_obs
    pct_vector = np.sign(adata_subs.X.toarray()).sum(axis=0) / n_obs
    return pd.Series(pct_vector, index=adata_subs.var_names, name=f"pct_{layer}")


def rank_degenes(de_results, adata, marker_df, logfc=1.5, pct=0.3):
    """
    Filter and rank DE genes
    """
    ranked_genes = dict()
    for layer in de_results.keys():
        marker_list = marker_df[marker_df["layer"] == layer].index.to_list()
        de_genes = de_results[layer].loc[marker_list]
        de_genes["pct"] = calculate_layerpct(adata, layer).loc[marker_list]
        filtered_markers = (
            de_genes[(de_genes["logFC"] > logfc) & (de_genes["pct"] > pct)]
            .sort_values(by="F", ascending=False)
            .copy()
        )
        ranked_genes[layer] = filtered_markers
    return ranked_genes
