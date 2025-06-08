# Preprocess Workflow

This directory contains the preprocessing notebooks used to prepare raw spatial and single-cell datasets for downstream analysis:

- **make_adata.ipynb**: Import raw data files and construct initial AnnData objects.
- **make_adata_filtered.ipynb**: Apply quality control filters to remove low-quality spots or cells.
- **make_pseudobulk.ipynb**: Aggregate expression data into pseudobulk profiles by sample and annotation groups.
