import gseapy
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm.notebook import tqdm


def create_rank_df(
    adata,
    group,
    key,
    gene_col="names",
    score_col="scores",
    lfc_col="logfoldchanges",
    method="score",
    pval_cutoff=None,
):
    """
    Creates a DataFrame with ranked genes based on the results of DE (differential expression) analysis
    """
    # read results for DE
    de_res = sc.get.rank_genes_groups_df(
        adata, group=group, key=key, pval_cutoff=pval_cutoff
    )
    de_res = de_res.set_index(gene_col)
    # choose a method to create ranked gene list based on DE results
    if method == "score":
        score_values = de_res[score_col]
    elif method == "lfc_product":
        score_values = de_res["logfoldchanges"].abs() * de_res[score_col]
    else:
        raise ValueError("No such method")
    # create df with ranked genes
    gene_rank_df = score_values.sort_values(ascending=False).to_frame()
    return gene_rank_df


def enrich_celltypes(
    celltypes,
    gene_sets,
    adata=None,
    key=None,
    gene_col="names",
    score_col="scores",
    lfc_col="logfoldchanges",
    method="score",
    pval_cutoff=None,
    gsea_kw=None,
):
    """
    Find enrichment of gene_sets`for ranked gene lists
    """
    # set gsea parameters
    gsea_standart_kw = {
        "threads": 4,
        "min_size": 0,
        "max_size": 1000,
        "permutation_num": 1000,  # reduce number to speed up testing
        "outdir": None,  # don't write to disk
        "seed": 4,
        "verbose": True,
    }

    if gsea_kw:
        gsea_standart_kw.update(gsea_kw)

    celltype_list = celltypes if adata else celltypes.keys()

    # results list
    res_df_list = list()
    # enrichment for each celltype
    for celltype in tqdm(celltype_list):
        # create a DataFrame with ranked genes
        if adata:
            rank_df = create_rank_df(
                adata,
                group=celltype,
                key=key,
                gene_col=gene_col,
                score_col=score_col,
                lfc_col=lfc_col,
                method=method,
                pval_cutoff=pval_cutoff,
            )
        else:
            rank_df = (
                celltypes[celltype]
                .set_index(gene_col)[score_col]
                .sort_values(ascending=False)
                .to_frame()
            )
        # enrich gene_sets based on ranked genes list
        gsea_res = gseapy.prerank(rnk=rank_df, gene_sets=gene_sets, **gsea_standart_kw)
        # save enrichment results
        res_df = gsea_res.res2d
        res_df["celltype"] = celltype
        res_df_list.append(res_df)
    # create DataFrame with results
    enrich_res = pd.concat(res_df_list, axis=0)
    return enrich_res


def process_enrichment_df(enrich_df):
    """
    Process the DataFrame with enrichment results
    """
    # convert FDR to -log10(FDR)
    enrich_df["-log10(FDR)"] = -np.log10(enrich_df["FDR q-val"].astype(float) + 1e-3)
    # Add direction column
    enrich_df["direction"] = enrich_df["ES"].map(
        lambda x: "enriched" if x > 0 else "depleted"
    )
    # add significance category column
    enrich_df["significant"] = enrich_df["FDR q-val"].map(
        lambda x: "FDR < 0.05" if x < 0.05 else "FDR >= 0.05"
    )


def enrichment_plot(
    enrich_res, figsize=(7, 9), dpi=100, xlabel="Enrichment score", **kw_scatterplot
):
    # create figure object
    fig = plt.figure(figsize=figsize, dpi=dpi)

    # create scatterplot
    ax = sns.scatterplot(data=enrich_res.reset_index(), **kw_scatterplot)

    # modify label parameters
    ax.tick_params(labelsize=10)
    ax.legend(fontsize=15)
    ax.set_xlabel(xlabel, fontsize=15)
    ax.set_ylabel("", fontsize=10)
    ax.grid(False)
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
