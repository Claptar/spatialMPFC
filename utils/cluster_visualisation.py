import numpy as np
import pandas as pd
from scipy.signal import medfilt2d, wiener, convolve2d
from tqdm.notebook import tqdm
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import ticker


def apply_filter(
    adata,
    value_columns,
    sample_column,
    method="median",
    size=3,
    key="_filtered",
    new_columns=None,
):
    """
    Apply a spatial filter to the specified value columns in an AnnData object.
    """
    # filter method
    if method == "median":
        filter_func = medfilt2d
        kwargs = {"kernel_size": 3}
    elif method == "wiener":
        filter_func = wiener
        kwargs = {"mysize": size}
    elif method == "mean":
        filter_matrix = np.ones((size, size)) / size**2
        filter_func = convolve2d
        kwargs = {"in2": filter_matrix, "mode": "same"}
    else:
        raise ValueError("Invalid filter method. Choose 'median' or 'wiener'")

    # rotation matrix
    A = 1 / np.sqrt(2) * np.array([[1, -1], [1, 1]])

    # get sample list
    sample_list = adata.obs[sample_column].unique().to_list()
    results = list()

    for sample in tqdm(sample_list):
        # sub-sample data
        obs_mask = adata.obs[sample_column] == sample
        raw_value = adata.obs.loc[obs_mask, value_columns].copy()
        coordinates = adata.obsm["spatial"][obs_mask.values]

        # result df template
        result_df = raw_value.copy()
        result_df.columns = result_df.columns + key

        # center data to perform rotation
        x_max, y_max = coordinates.max(axis=0)
        x_min, y_min = coordinates.min(axis=0)

        x_center = x_min + (x_max - x_min) / 2
        y_center = y_min + (y_max - y_min) / 2

        centered_coord = coordinates - np.array([x_center, y_center])

        # rotate data
        rotate_coord = centered_coord @ A

        # convert coordinates to grid coordinates
        floor_divide = np.floor_divide(rotate_coord, 200)
        x_max, y_max = floor_divide.max(axis=0)
        x_min, y_min = floor_divide.min(axis=0)
        real_coord = floor_divide - np.array([x_min, y_min])

        # create 2d matrix template for the data
        x = np.arange(0, x_max - x_min + 1)
        y = np.arange(0, y_max - y_min + 1)
        coord = np.array(np.meshgrid(x, y)).T.reshape(-1, 2)

        # create 2d representation of the data
        df = pd.DataFrame(coord, columns=["x", "y"]).set_index(["x", "y"])
        df[value_columns] = 0
        df.loc[real_coord.tolist(), value_columns] = raw_value.values

        for value in df.columns.to_list():
            X = df.reset_index().pivot(index="x", columns="y", values=value).values
            X_filt = filter_func(X, **kwargs)
            result_df[value + key] = (
                pd.DataFrame(X_filt).T.unstack().loc[real_coord.tolist()].values
            )
        results.append(result_df)

    # concat data
    filtered_df = pd.concat(results, axis=0)
    new_columns = filtered_df.columns if new_columns is None else new_columns
    adata.obs[new_columns] = filtered_df


def plot_cluster_heatmap(
    df,
    order,
    labels,
    condition_series,
    clusters,
    cluster_cmap="Set3",
    cond_cmap="tab10",
    cond_palette=None,
    figsize=(10, 8),
    vmax=0.3,
    vmin=-0.3,
    center=0,
    ax=None,
    fig=None,
    cbar=True,
):
    """
    Plot a heatmap of df[order] with top cluster bar and left condition bar.

    Parameters:
    - df: DataFrame (samples x features)
    - order: list of features for column ordering
    - labels: Series mapping feature -> cluster label
    - condition_series: pandas Series (sample index -> condition label)
    - clusters: iterable of cluster labels in desired order
    - cluster_cmap: colormap name for clusters
    - cond_cmap: colormap name for conditions
    - figsize: tuple for figure size
    - vmax, vmin, center: heatmap scale parameters

    Returns:
    - fig, ax: matplotlib Figure and Axes
    """
    mpl.rcdefaults()
    if ax is None or fig is None:
        fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(
        df[order],
        ax=ax,
        cmap="RdBu_r",
        vmax=vmax,
        vmin=vmin,
        center=center,
        cbar=cbar,
    )
    ax.set_yticks([])
    divider = make_axes_locatable(ax)

    # Top cluster bar
    cax_top = divider.append_axes("top", size="2%", pad=0.05)
    cmap1 = plt.get_cmap(cluster_cmap)
    sizes1 = labels.value_counts().loc[clusters].values
    pos1 = sizes1.cumsum()
    bounds1 = [0] + pos1.tolist()
    norm1 = mpl.colors.BoundaryNorm(bounds1, cmap1.N)
    fig.colorbar(
        mpl.cm.ScalarMappable(cmap=cmap1, norm=norm1),
        cax=cax_top,
        ticks=bounds1,
        orientation="horizontal",
        spacing="proportional",
    )
    cax_top.xaxis.set_major_locator(ticker.FixedLocator(pos1 - sizes1 / 2))
    cax_top.xaxis.set_major_formatter(ticker.FixedFormatter(clusters))
    cax_top.xaxis.tick_top()

    # Left condition bar
    cax_left = divider.append_axes("left", size="2%", pad=0.05)
    cmap2 = (
        plt.get_cmap(cond_cmap)
        if cond_palette is None
        else mpl.colors.ListedColormap(list(cond_palette.values())[::-1])
    )
    cond_counts = (
        condition_series.value_counts()
        if cond_palette is None
        else condition_series.value_counts()
    )
    conds = cond_counts.index[::-1]
    sizes2 = cond_counts.loc[conds].values
    pos2 = sizes2.cumsum()
    bounds2 = [0] + pos2.tolist()
    norm2 = mpl.colors.BoundaryNorm(bounds2, cmap2.N)
    fig.colorbar(
        mpl.cm.ScalarMappable(cmap=cmap2, norm=norm2),
        cax=cax_left,
        ticks=bounds2,
        orientation="vertical",
        spacing="proportional",
    )
    cax_left.yaxis.set_major_locator(ticker.FixedLocator(pos2 - sizes2 / 3))
    cax_left.yaxis.set_major_formatter(ticker.FixedFormatter(conds))
    cax_left.yaxis.tick_left()
    cax_left.set_yticklabels(cax_left.get_yticklabels(), rotation=90, fontsize=12)

    return fig, ax


def plot_cluster_profiles(
    layer_mean,
    labels,
    layers,
    specie_palette,
    n_clusters=None,
    figsize=(14, None),
    hspace=0.4,
    marker=".",
    title_prefix="Cluster ",
):
    """
    Plot average expression profiles for each cluster across layers and conditions.

    Parameters:
    - layer_mean: DataFrame with MultiIndex (condition, layer) and genes as columns
    - labels: Series mapping gene names to cluster labels
    - layers: list of layer names in order
    - specie_palette: dict mapping condition names to colors
    - n_clusters: number of clusters; defaults to labels.nunique()
    - figsize: (width, height) tuple or (width, None) to auto-calc height per cluster
    - hspace: vertical spacing between subplots
    - marker: marker style for line plots
    - title_prefix: prefix for each subplot title

    Returns:
    - fig, axes: matplotlib figure and array of Axes
    """
    mpl.rcdefaults()
    if n_clusters is None:
        n_clusters = labels.nunique()
    nrows = (n_clusters + 1) // 2
    # auto height: 3.5 per row if not specified
    height = figsize[1] or (3.5 * nrows)
    fig, axes = plt.subplots(
        nrows,
        2,
        figsize=(figsize[0], height),
        gridspec_kw={"hspace": hspace},
    )

    for label, ax in zip(sorted(labels.unique()), axes.flatten()):
        genes = labels[labels == label].index.tolist()
        # mean across genes per (condition, layer)
        cluster_mean = (
            layer_mean[genes]
            .mean(axis=1)
            .reorder_levels(["layer", "condition"])
            .unstack()
            .loc[layers]
        )
        cluster_mean.plot.line(color=specie_palette, ax=ax, marker=marker)

        ax.legend(
            fontsize=10,
            loc="upper left",
            bbox_to_anchor=(1, 1),
            prop={"size": 12},
        )
        if label % 2 == 1:
            ax.get_legend().remove()
        ax.set_title(f"{title_prefix}{label}", fontsize=16)
        ax.grid(False)

    # hide empty subplot if odd number
    if n_clusters % 2 == 1:
        axes.flatten()[-1].axis("off")

    return fig, axes


def plot_cluster_spline_profiles(
    grid,
    df_spline,
    labels,
    layers,
    specie_palette,
    n_clusters=None,
    figsize=(14, None),
    hspace=0.4,
    title_prefix="Cluster ",
):
    """
    Plot spline-approximated expression profiles for each cluster across layers and conditions.

    Parameters:
    - grid: array-like x-axis values (e.g., continuous layer grid)
    - df_spline: DataFrame containing spline values with columns for each gene and a 'condition' column
    - labels: Series mapping gene names to cluster labels
    - layers: list of discrete layer labels for x-ticks
    - specie_palette: dict mapping condition names to colors
    - n_clusters: number of clusters (defaults to labels.nunique())
    - figsize: (width, height) tuple or (width, None) to auto-calc height
    - hspace: vertical space between rows
    - title_prefix: prefix for subplot titles

    Returns:
    - fig, axes: matplotlib Figure and Axes array
    """

    mpl.rcdefaults()
    if n_clusters is None:
        n_clusters = labels.nunique()
    nrows = (n_clusters + 1) // 2
    height = figsize[1] or (3.5 * nrows)
    fig, axes = plt.subplots(
        nrows, 2, figsize=(figsize[0], height), gridspec_kw={"hspace": hspace}
    )

    for label, ax in zip(sorted(labels.unique()), axes.flatten()):
        genes = labels[labels == label].index.tolist()
        # compute cluster mean spline per condition
        df_cluster = df_spline[genes].mean(axis=1).to_frame(name="value")
        df_cluster["condition"] = df_spline["condition"]
        df_cluster_pivot = df_cluster.pivot(columns="condition", values="value")
        df_cluster_pivot.index = grid
        df_cluster_pivot.plot.line(color=specie_palette, ax=ax)

        ax.legend(
            fontsize=10, loc="upper left", bbox_to_anchor=(1, 1), prop={"size": 12}
        )
        ax.xaxis.set_major_locator(ticker.FixedLocator(np.arange(1, len(layers) + 1)))
        ax.xaxis.set_major_formatter(ticker.FixedFormatter(layers))
        if label % 2 == 1:
            ax.get_legend().remove()
        ax.set_title(f"{title_prefix}{label}", fontsize=16)
        ax.grid(False)

    # hide empty subplot if odd count
    if n_clusters % 2 == 1:
        axes.flatten()[-1].axis("off")

    return fig, axes
