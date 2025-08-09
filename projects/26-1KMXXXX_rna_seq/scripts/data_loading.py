
import os
import scanpy as sc
import pandas as pd
import numpy as np

# Data Loading
path_GSM8440146 = 'data/GSE273941_RAW/GSM8440146'
path_GSM8440147 = 'data/GSE273941_RAW/GSM8440147'

#print("Files in GSM8440146 (control):")
#print(os.listdir(path_GSM8440146))

#print("Files in GSM8440147 (ACM treated):")
#print(os.listdir(path_GSM8440147))

adata_control = sc.read_10x_mtx(path_GSM8440146, cache=True)
adata_treated = sc.read_10x_mtx(path_GSM8440147, cache=True)


#print(adata_control)
#print(adata_treated)

# 1. Annotate sample origin
adata_control.obs["condition"] = "control"
adata_treated.obs["condition"] = "treated"

# 2. Merge datasets
adata = adata_control.concatenate(adata_treated, batch_key="batch", batch_categories=["control", "treated"])

# Filter cells with too few genes or high mitochondrial content
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

adata = adata[adata.obs['n_genes_by_counts'] > 200, :]
adata = adata[adata.obs['pct_counts_mt'] < 10, :]

# Group cells by sample
pseudobulk_counts = (
    adata.to_df()
    .groupby(adata.obs['sample'])  # this will be 'control_1', 'treated_1'
    .sum()
    .T  # Transpose → genes × samples
)

# Save for DESeq2
pseudobulk_counts.to_csv("pseudobulk_counts.csv")

meta = pd.DataFrame({
    'sample': ['control_1', 'treated_1'],
    'condition': ['control', 'treated']
})
meta.to_csv("sample_metadata.csv", index=False)

""""
# 3. Basic QC metrics
adata.var_names_make_unique()
adata.obs["n_counts"] = adata.X.sum(axis=1).A1  # Total counts per cell
adata.obs["n_genes"] = (adata.X > 0).sum(axis=1).A1  # Number of expressed genes per cell

# 4. Calculate % mitochondrial genes (optional but often useful)
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

# 5. Filter cells (low quality / dead / doublets)
sc.pp.filter_cells(adata, min_genes=200)
adata = adata[adata.obs["pct_counts_mt"] < 10, :]

# 6. Filter genes (remove rarely expressed genes)
sc.pp.filter_genes(adata, min_cells=3)

# 7. Library size normalization (for model offset)
# Compute size factors (library size per cell)
adata.obs["size_factors"] = adata.obs["total_counts"] / np.median(adata.obs["total_counts"])

# NOTE: Don't log1p-transform here – we want raw counts for modeling
# But you can keep a copy for PCA etc. if needed:
adata.raw = adata.copy()

# 8. Optional normalization/log1p for visualization (not for modeling)
# sc.pp.normalize_total(adata, target_sum=1e4)
# sc.pp.log1p(adata)

# 9. Save intermediate state
adata.write("adata_preprocessed.h5ad")

# 10. Final check
print(adata)
print(adata.obs[["condition", "n_counts", "n_genes", "pct_counts_mt", "size_factors"]].head())
print(adata.obs["condition"])
"""

"""
# Data Exploration

adata_control.var_names_make_unique()
adata_control.var['mt'] = adata_control.var_names.str.upper().str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata_control, qc_vars=['mt'], inplace=True)
print(adata_control.obs.head())

import matplotlib.pyplot as plt
sc.pl.violin(adata_control, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4)

sc.pp.filter_genes(adata_control, min_cells=3)

adata_control = adata_control[adata_control.obs['n_genes_by_counts'] > 200, :]
adata_control = adata_control[adata_control.obs['pct_counts_mt'] < 5, :]

sc.pp.normalize_total(adata_control, target_sum=1e4)
sc.pp.log1p(adata_control)
"""