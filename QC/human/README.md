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
os.listdir('../../data/anndata_objects/human')
```




    ['human.h5ad',
     'human_759.h5ad',
     'human_j12.h5ad',
     'human_j3.h5ad',
     'human_j4.h5ad',
     'human_j6.h5ad']




```python
adata = sc.read_h5ad('../../data/anndata_objects/human/human.h5ad')
adata.obs_names_make_unique()
adata
```

    C:\Users\aleks\anaconda3\envs\scanorama39\lib\site-packages\anndata\_core\anndata.py:1828: UserWarning: Observation names are not unique. To make them unique, call `.obs_names_make_unique`.
      utils.warn_names_duplicates("obs")
    




    AnnData object with n_obs × n_vars = 17633 × 19966
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




    WM             5004
    L3             3118
    L5             1951
    L4             1869
    Empty spots    1411
    L1             1389
    L6             1167
    L2             1075
    L6b             339
    L6a             309
    Name: label, dtype: int64




```python
samples = list(adata.uns['spatial'].keys())
samples
```




    ['human_759', 'human_j12', 'human_j3', 'human_j4', 'human_j6']



## Посмотрим на изображения


```python
sq.pl.spatial_scatter(adata, library_key='sample_id', title=samples, ncols=5)
```


    
![png](./figures/output_10_0.png)
    



```python
sq.pl.spatial_scatter(adata, color=['label'], library_key='sample_id', size=1.1, ncols=5, title=samples)
```


    
![png](./figures/output_11_0.png)
    


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
        finished (0:02:14)
    


```python
sc.pl.pca(adata, 
          color=['label', 'sample_id'],
          ncols=2, show=True)
```

    C:\Users\aleks\anaconda3\envs\scanorama39\lib\site-packages\scanpy\plotting\_tools\scatterplots.py:392: UserWarning: No data for colormapping provided via 'c'. Parameters 'cmap' will be ignored
      cax = scatter(
    C:\Users\aleks\anaconda3\envs\scanorama39\lib\site-packages\scanpy\plotting\_tools\scatterplots.py:392: UserWarning: No data for colormapping provided via 'c'. Parameters 'cmap' will be ignored
      cax = scatter(




<img src="./figures/output_19_1.png" width="70%" height="70%" />
    


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

    C:\Users\aleks\AppData\Local\Temp\ipykernel_24060\737622212.py:1: ImplicitModificationWarning: Trying to modify attribute `.obs` of view, initializing view as actual.
      adata.obs['cluster'] = adata.obs.label
    




    WM             5004
    L3             3118
    L5             1951
    L4             1869
    L6             1815
    Empty spots    1411
    L1             1389
    L2             1075
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
        'pvals_adj', sorted np.recarray to be indexed by group ids (0:01:51)
    


```python
sc.pl.rank_genes_groups_dotplot(adata[adata.obs.cluster != 'Empty spots'],
                                groupby="cluster", n_genes=5, values_to_plot='logfoldchanges',
                                min_logfoldchange=1, vmax=3, vmin=-3, cmap='bwr', save='dotplot_fc_human.pdf')
```

    WARNING: dendrogram data not found (using key=dendrogram_cluster). Running `sc.tl.dendrogram` with default parameters. For fine tuning it is recommended to run `sc.tl.dendrogram` independently.
        using 'X_pca' with n_pcs = 50
    Storing dendrogram info using `.uns['dendrogram_cluster']`
    

    C:\Users\aleks\anaconda3\envs\scanorama39\lib\site-packages\anndata\compat\_overloaded_dict.py:106: ImplicitModificationWarning: Trying to modify attribute `._uns` of view, initializing view as actual.
      self.data[key] = value
    

    WARNING: saving figure to file figures\dotplot_dotplot_fc_human.pdf
    

    C:\Users\aleks\anaconda3\envs\scanorama39\lib\site-packages\scanpy\plotting\_dotplot.py:749: UserWarning: No data for colormapping provided via 'c'. Parameters 'cmap', 'norm' will be ignored
      dot_ax.scatter(x, y, **kwds)
    


    
![png](./figures/output_24_4.png)
    



```python
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, groupby="cluster", min_logfoldchange=1, save='dotplot_expr_human.pdf')
```

    WARNING: dendrogram data not found (using key=dendrogram_cluster). Running `sc.tl.dendrogram` with default parameters. For fine tuning it is recommended to run `sc.tl.dendrogram` independently.
        using 'X_pca' with n_pcs = 50
    Storing dendrogram info using `.uns['dendrogram_cluster']`
    WARNING: Groups are not reordered because the `groupby` categories and the `var_group_labels` are different.
    categories: Empty spots, L1, L2, etc.
    var_group_labels: L1, WM, L6, etc.
    WARNING: saving figure to file figures\dotplot_dotplot_expr_human.pdf
    

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
    Finish (0:00:03)
    


```python
sq.gr.nhood_enrichment(adata, cluster_key="cluster")
sq.pl.nhood_enrichment(adata, cluster_key="cluster", method="average")
```

    Calculating neighborhood enrichment using `1` core(s)
    


      0%|          | 0/1000 [00:00<?, ?/s]


    Adding `adata.uns['cluster_nhood_enrichment']`
    Finish (0:00:15)




<img src="./figures/output_28_3.png" width="50%" height="50%" />
    



```python
sq.gr.spatial_autocorr(adata, mode="moran", genes=adata.var_names)
```

    C:\Users\aleks\anaconda3\envs\scanorama39\lib\site-packages\scanpy\metrics\_gearys_c.py:293: UserWarning: 2233 variables were constant, will return nan for these.
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
      <th>MTRNR2L12</th>
      <td>0.743979</td>
      <td>0.0</td>
      <td>0.000019</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>MTRNR2L8</th>
      <td>0.728145</td>
      <td>0.0</td>
      <td>0.000019</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>NRGN</th>
      <td>0.578565</td>
      <td>0.0</td>
      <td>0.000019</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>PTGDS</th>
      <td>0.537307</td>
      <td>0.0</td>
      <td>0.000019</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>TUBA1A</th>
      <td>0.524920</td>
      <td>0.0</td>
      <td>0.000019</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>MT3</th>
      <td>0.522044</td>
      <td>0.0</td>
      <td>0.000019</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>SPARCL1</th>
      <td>0.498617</td>
      <td>0.0</td>
      <td>0.000019</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>COL1A1</th>
      <td>0.485158</td>
      <td>0.0</td>
      <td>0.000019</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>NEFL</th>
      <td>0.483981</td>
      <td>0.0</td>
      <td>0.000019</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>COL1A2</th>
      <td>0.467839</td>
      <td>0.0</td>
      <td>0.000019</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>COL3A1</th>
      <td>0.467472</td>
      <td>0.0</td>
      <td>0.000019</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>CST3</th>
      <td>0.453844</td>
      <td>0.0</td>
      <td>0.000019</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>HSP90AA1</th>
      <td>0.453142</td>
      <td>0.0</td>
      <td>0.000019</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>HBG2</th>
      <td>0.446280</td>
      <td>0.0</td>
      <td>0.000019</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>APOD</th>
      <td>0.438664</td>
      <td>0.0</td>
      <td>0.000019</td>
      <td>NaN</td>
    </tr>
  </tbody>
</table>
</div>




```python
sq.pl.spatial_scatter(adata,
                      color=['MTRNR2L12', 'CXCL14', 'ENC1', 'label'],
                      library_key='sample_id', ncols=4, img=True, size=1.3,
                      save='../figures/markers_human.pdf')
```


    
![png](./figures/output_31_0.png)
    

