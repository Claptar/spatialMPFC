```python
import warnings
import scanpy as sc
import squidpy as sq
import anndata as an
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from urllib import request
import scanorama
import json
import os

sc.settings.set_figure_params(dpi=80)
#sc.set_figure_params(facecolor="white", figsize=(8, 8))
warnings.simplefilter(action='ignore', category=FutureWarning)
sc.settings.verbosity = 3
```

# Глобальное окружение

## Подгрузим данные


```python
os.listdir('../../data/anndata_objects/chimp')
```




    ['chimp.h5ad',
     'chimp_11454.h5ad',
     'chimp_13302.h5ad',
     'chimp_13309.h5ad',
     'chimp_j11.h5ad',
     'chimp_j8.h5ad']




```python
adata = sc.read_h5ad('../../data/anndata_objects/chimp/chimp.h5ad')
adata.obs_names_make_unique()
adata
```

    C:\Users\aleks\anaconda3\envs\scanorama39\lib\site-packages\anndata\_core\anndata.py:1828: UserWarning: Observation names are not unique. To make them unique, call `.obs_names_make_unique`.
      utils.warn_names_duplicates("obs")
    




    AnnData object with n_obs × n_vars = 22088 × 25080
        obs: 'in_tissue', 'array_row', 'array_col', 'label', 'sample_id'
        var: 'gene_ids', 'feature_types'
        uns: 'spatial'
        obsm: 'spatial'



Оставим споты только с разметкой


```python
adata = adata[adata.obs.label.notna()]
```


```python
adata.obs.label.value_counts()
```




    WM             6397
    L6             4710
    L3             3400
    L5             2870
    Empty spots    1726
    L4             1429
    L1             1024
    L2              530
    Name: label, dtype: int64




```python
samples = list(adata.uns['spatial'].keys())
samples
```




    ['chimp_11454', 'chimp_13302', 'chimp_13309', 'chimp_j11', 'chimp_j8']



## Посмотрим на изображения


```python
sq.pl.spatial_scatter(adata, library_key='sample_id', title=samples, ncols=5)
```


    
![png](./figures/output_10_0.png)
    



```python
sq.pl.spatial_scatter(adata, color=['label'], library_key='sample_id', size=1.1, ncols=5, title=samples)
```

    C:\Users\aleks\anaconda3\envs\scanorama39\lib\site-packages\anndata\compat\_overloaded_dict.py:106: ImplicitModificationWarning: Trying to modify attribute `._uns` of view, initializing view as actual.
      self.data[key] = value
    


    
![png](./figures/output_11_1.png)
    


## Посмотрим на метрики


```python
# add info on mitochondrial and hemoglobin genes to the objects.
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['hb'] = adata.var_names.str.contains(("^HB[AB]"))
adata.var['ribo'] = adata.var_names.str.contains(("^RP[LS]"))

sc.pp.calculate_qc_metrics(
    adata, qc_vars=['mt','hb','ribo'],
    percent_top=None, log1p=False, inplace=True)
```


```python
sc.pl.violin(adata,
             keys = ['n_genes_by_counts', 'total_counts','pct_counts_ribo', 'pct_counts_mt'],
             jitter=0.4, groupby = 'sample_id', rotation= 45)
```


    
![png](./figures/output_14_0.png)
    



```python
sc.pl.violin(adata,
             keys = ['n_genes_by_counts', 'total_counts','pct_counts_ribo', 'pct_counts_mt'],
             jitter=0.4, groupby = 'label', rotation= 90)
```


    
![png](./figures/output_15_0.png)
    



```python
sq.pl.spatial_scatter(adata,
                      color=['n_genes_by_counts', 'total_counts','pct_counts_ribo', 'pct_counts_mt', 'label'],
                      library_key='sample_id', ncols=5, img=False, size=1.3)
```


    
![png](./figures/output_16_0.png)
    


## PCA


```python
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.tl.pca(adata)
```

    normalizing counts per cell
        finished (0:00:00)
    computing PCA
        with n_comps=50
        finished (0:02:29)
    


```python
sc.pl.pca(adata, 
          color=['label', 'sample_id'],
          ncols=2, show=True)
```

    C:\Users\aleks\anaconda3\envs\scanorama39\lib\site-packages\scanpy\plotting\_tools\scatterplots.py:392: UserWarning: No data for colormapping provided via 'c'. Parameters 'cmap' will be ignored
      cax = scatter(
    C:\Users\aleks\anaconda3\envs\scanorama39\lib\site-packages\scanpy\plotting\_tools\scatterplots.py:392: UserWarning: No data for colormapping provided via 'c'. Parameters 'cmap' will be ignored
      cax = scatter(
    


    
![png](./figures/output_19_1.png)
    


### Дифф. экспрессия


```python
remove = adata.var_names.str.startswith('MT-')
keep = np.invert(remove)
print(sum(remove))

adata = adata[:,keep]
```

    13
    


```python
adata.obs['cluster'] = adata.obs.label
adata.obs.cluster.replace({'L6a': 'L6', 'L6b': 'L6'}, inplace=True)
adata.obs.cluster.value_counts()
```

    C:\Users\aleks\AppData\Local\Temp\ipykernel_11312\737622212.py:1: ImplicitModificationWarning: Trying to modify attribute `.obs` of view, initializing view as actual.
      adata.obs['cluster'] = adata.obs.label
    




    WM             6397
    L6             4710
    L3             3400
    L5             2870
    Empty spots    1726
    L4             1429
    L1             1024
    L2              530
    Name: cluster, dtype: int64




```python
sc.tl.rank_genes_groups(adata, groupby='cluster', method='wilcoxon', groups=['L1', 'WM', 'L6', 'L2', 'L3', 'L4', 'L5'])
```

    ranking genes
        finished: added to `.uns['rank_genes_groups']`
        'names', sorted np.recarray to be indexed by group ids
        'scores', sorted np.recarray to be indexed by group ids
        'logfoldchanges', sorted np.recarray to be indexed by group ids
        'pvals', sorted np.recarray to be indexed by group ids
        'pvals_adj', sorted np.recarray to be indexed by group ids (0:03:12)
    


```python
sc.pl.rank_genes_groups_dotplot(adata[adata.obs.cluster != 'Empty spots'],
                                groupby="cluster", n_genes=5, values_to_plot='logfoldchanges',
                                min_logfoldchange=1, vmax=3, vmin=-3, cmap='bwr')
```

    WARNING: dendrogram data not found (using key=dendrogram_cluster). Running `sc.tl.dendrogram` with default parameters. For fine tuning it is recommended to run `sc.tl.dendrogram` independently.
        using 'X_pca' with n_pcs = 50
    Storing dendrogram info using `.uns['dendrogram_cluster']`
    

    C:\Users\aleks\anaconda3\envs\scanorama39\lib\site-packages\anndata\compat\_overloaded_dict.py:106: ImplicitModificationWarning: Trying to modify attribute `._uns` of view, initializing view as actual.
      self.data[key] = value
    C:\Users\aleks\anaconda3\envs\scanorama39\lib\site-packages\scanpy\plotting\_dotplot.py:749: UserWarning: No data for colormapping provided via 'c'. Parameters 'cmap', 'norm' will be ignored
      dot_ax.scatter(x, y, **kwds)
    

    WARNING: saving figure to file figures\dotplot_dotplot_fc_human.pdf
    


    
![png](./figures/output_24_3.png)
    



```python
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, groupby="cluster", min_logfoldchange=1)
```

    WARNING: dendrogram data not found (using key=dendrogram_cluster). Running `sc.tl.dendrogram` with default parameters. For fine tuning it is recommended to run `sc.tl.dendrogram` independently.
        using 'X_pca' with n_pcs = 50
    Storing dendrogram info using `.uns['dendrogram_cluster']`
    WARNING: Groups are not reordered because the `groupby` categories and the `var_group_labels` are different.
    categories: Empty spots, L1, L2, etc.
    var_group_labels: L1, WM, L6, etc.
    

    C:\Users\aleks\anaconda3\envs\scanorama39\lib\site-packages\scanpy\plotting\_dotplot.py:749: UserWarning: No data for colormapping provided via 'c'. Parameters 'cmap', 'norm' will be ignored
      dot_ax.scatter(x, y, **kwds)
    


    
![png](./figures/output_25_2.png)
    


## Spatial features


```python
sq.gr.spatial_neighbors(
    adata, coord_type="generic", library_key="sample_id", delaunay=True
)
```

    Creating graph using `generic` coordinates and `None` transform and `5` libraries.
    Adding `adata.obsp['spatial_connectivities']`
           `adata.obsp['spatial_distances']`
           `adata.uns['spatial_neighbors']`
    Finish (0:00:09)
    


```python
sq.gr.nhood_enrichment(adata, cluster_key="cluster")
sq.pl.nhood_enrichment(adata, cluster_key="cluster", method="average")
```

    Calculating neighborhood enrichment using `1` core(s)
    


      0%|          | 0/1000 [00:00<?, ?/s]


    Adding `adata.uns['cluster_nhood_enrichment']`
    Finish (0:00:24)
    


    
![png](./figures/output_28_3.png)
    



```python
sq.gr.spatial_autocorr(adata, mode="moran", genes=adata.var_names)
```

    C:\Users\aleks\anaconda3\envs\scanorama39\lib\site-packages\scanpy\metrics\_gearys_c.py:293: UserWarning: 6598 variables were constant, will return nan for these.
      warnings.warn(
    

    Calculating moran's statistic for `None` permutations using `1` core(s)
    Adding `adata.uns['moranI']`
    Finish (0:00:00)
    


```python
adata.uns["moranI"].head(15)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>I</th>
      <th>pval_norm</th>
      <th>var_norm</th>
      <th>pval_norm_fdr_bh</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>MBP</th>
      <td>0.620382</td>
      <td>0.0</td>
      <td>0.000015</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>ENSPTRG00000047423</th>
      <td>0.564648</td>
      <td>0.0</td>
      <td>0.000015</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>PLP1</th>
      <td>0.554664</td>
      <td>0.0</td>
      <td>0.000015</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>GFAP</th>
      <td>0.533049</td>
      <td>0.0</td>
      <td>0.000015</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>HBA1</th>
      <td>0.497536</td>
      <td>0.0</td>
      <td>0.000015</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>HBB</th>
      <td>0.422591</td>
      <td>0.0</td>
      <td>0.000015</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>HBG2</th>
      <td>0.415482</td>
      <td>0.0</td>
      <td>0.000015</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>NPY</th>
      <td>0.408818</td>
      <td>0.0</td>
      <td>0.000015</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>HBA2</th>
      <td>0.370459</td>
      <td>0.0</td>
      <td>0.000015</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>THBS1</th>
      <td>0.368609</td>
      <td>0.0</td>
      <td>0.000015</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>IGF2</th>
      <td>0.360383</td>
      <td>0.0</td>
      <td>0.000015</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>HBG1</th>
      <td>0.355704</td>
      <td>0.0</td>
      <td>0.000015</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>TUBA1A</th>
      <td>0.350845</td>
      <td>0.0</td>
      <td>0.000015</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>NEFL</th>
      <td>0.350207</td>
      <td>0.0</td>
      <td>0.000015</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>CALB2</th>
      <td>0.321590</td>
      <td>0.0</td>
      <td>0.000015</td>
      <td>NaN</td>
    </tr>
  </tbody>
</table>
</div>




```python
sq.pl.spatial_scatter(adata,
                      color=['VIM', 'GFAP', 'CCK', 'label'],
                      library_key='sample_id', ncols=4, img=True, size=1.3,
                      save='../figures/markers_human.pdf')
```


    
![png](./figures/output_31_0.png)
    

