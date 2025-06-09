import numpy as np
import requests
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from pandas.api.types import CategoricalDtype
from matplotlib import rcdefaults


def p_val_group(value):
    """
    Classify p-value into significance categories.
    """
    if value >= 0.05:
        return "No significance"
    elif 0.01 <= value < 0.05:
        return "p < 0.05"
    elif 0.001 <= value < 0.01:
        return "p < 0.01"
    else:
        return "p < 0.001"


class GeneCluster:
    """
    Class representing a gene cluster for functional enrichment analysis.
    """

    def __init__(self, genes, descr, label, background_genes, **kwargs):
        self.label = label
        self.genes = genes
        self.descr = descr
        self.background_genes = background_genes
        self.userlist_id = self._get_userlist_id()
        self.background_id = self._get_background_id()
        self.enrichment_res = dict()

    def _get_userlist_id(self):
        """
        Create a user list on SpeedRichR with the provided genes and description.
        """
        base_url = "https://maayanlab.cloud/speedrichr"

        description = "sample gene set with background"

        res = requests.post(
            base_url + "/api/addList",
            files=dict(
                list=(None, "\n".join(self.genes)),
                description=(None, description),
            ),
        )
        if res.ok:
            userlist_response = res.json()
        else:
            raise Exception("Error analyzing gene list")
        return userlist_response["userListId"]

    def _get_background_id(self):
        """
        Create a background list on SpeedRichR with the provided background genes.
        """
        base_url = "https://maayanlab.cloud/speedrichr"

        res = requests.post(
            base_url + "/api/addbackground",
            data=dict(background="\n".join(self.background_genes)),
        )

        if res.ok:
            background_response = res.json()
        else:
            raise Exception("Error analyzing gene list")
        return background_response["backgroundid"]

    def enrich(self, gene_set_library):
        """
        Perform functional enrichment analysis using SpeedRichR for the specified gene set library.
        """
        # get enrichment data/human_specific_genes
        base_url = "https://maayanlab.cloud/speedrichr"

        res = requests.post(
            base_url + "/api/backgroundenrich",
            data=dict(
                userListId=self.userlist_id,
                backgroundid=self.background_id,
                backgroundType=gene_set_library,
            ),
        )
        if res.ok:
            data = res.json()
        else:
            raise Exception("Error analyzing gene list")
        # convert data/human_specific_genes to df
        columns = [
            "Rank",
            "Term",
            "p-val",
            "Z-score",
            "Combined score",
            "Overlapping genes",
            "Adjusted P-value",
            "Old p-value",
            "Old adjusted p-value",
        ]
        results = pd.DataFrame(data[gene_set_library], columns=columns)
        # preprocess df
        results.Term = results.Term.astype(str)
        results["num_overlap_genes"] = results["Overlapping genes"].apply(
            lambda x: len(x)
        )
        results["neg_log10(p_adj)"] = -np.log10(results["Adjusted P-value"])
        results["cluster_label"] = self.label
        # save to enrichment_res
        self.enrichment_res[gene_set_library] = results

    @staticmethod
    def enrich_geneclusters(geneclusters, gene_set_library):
        """
        Perform functional enrichment analysis for a list of GeneCluster objects.
        """
        for gc in geneclusters:
            gc.enrich(gene_set_library)


def scatter_enrichment(
    enrich_res: pd.DataFrame,
    terms: list,
    color_discrete_map: dict,
    db_name: str,
    *,
    pval_groups=("No significance", "p < 0.05", "p < 0.01", "p < 0.001"),
    figsize=(7, 13),
    dpi=100,
    size_range=(20, 250),
    fontsize_row=5,
    fontsize_axis=15,
):
    """
    Draw a p-value coloured scatter plot of enrichment results.

    Parameters
    ----------
    enrich_res : DataFrame
        Full enrichment results table (must contain the columns
        'Adjusted P-value', 'cluster_label', 'Term', 'num_overlap_genes').
    terms : list[str]
        Sub-set of term IDs/indices to plot (row labels of `enrich_res`).
    utils : module
        Must expose `p_val_group(p)`  -> str, returning one of `pval_groups`.
    color_discrete_map : dict
        Mapping from p-value group labels to colours.
    db_name : str
        Title that will appear on the plot.
    pval_groups : tuple[str], optional
        Allowed p-value group labels in ascending significance.
    figsize, dpi, size_range, fontsize_row, fontsize_axis : misc
        Aesthetics.

    Returns
    -------
    fig, ax : (Figure, Axes)
    """
    # ---- 1. filter & annotate ------------------------------------------------
    enrich_res_plot = enrich_res.loc[terms].copy()

    cat_type = CategoricalDtype(categories=pval_groups, ordered=True)
    enrich_res_plot["p-value"] = (
        enrich_res_plot["Adjusted P-value"].apply(p_val_group).astype(cat_type)
    )

    # ---- 2. create the figure ------------------------------------------------
    rcdefaults()  # reset to Matplotlib defaults
    fig = plt.figure(figsize=figsize, dpi=dpi)

    ax = sns.scatterplot(
        data=enrich_res_plot.reset_index(),
        x="cluster_label",
        y="Term",
        size="num_overlap_genes",
        sizes=size_range,
        hue="p-value",
        palette=color_discrete_map,
    )

    # ---- 3. styling ----------------------------------------------------------
    ax.tick_params(labelsize=10)
    ax.legend(fontsize=fontsize_axis)
    ax.set_xlabel("Cluster label", fontsize=fontsize_axis)
    ax.set_ylabel("", fontsize=fontsize_row)
    ax.grid(False)
    ax.set_title(db_name)
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))

    return fig, ax
