import os
import scanpy as sc
import pandas as pd
import numpy as np

# === 1. Data Loading ===
path_GSM8440146 = 'data/GSE273941_RAW/GSM8440146'  # Control
path_GSM8440147 = 'data/GSE273941_RAW/GSM8440147'  # Treated

adata_control = sc.read_10x_mtx(path_GSM8440146, cache=True)
adata_treated = sc.read_10x_mtx(path_GSM8440147, cache=True)

# === 2. Annotate sample origin ===
adata_control.obs["condition"] = "control"
adata_treated.obs["condition"] = "treated"

# === 3. Merge datasets ===
adata = adata_control.concatenate(
    adata_treated,
    batch_key="batch",
    batch_categories=["control", "treated"]
)

# Optional: Assign simplified sample IDs for pseudobulk
adata.obs["sample"] = adata.obs["batch"].map({
    "control": "control_1",
    "treated": "treated_1"
})

# === 4. QC filtering ===
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

# Basic filters (can be tuned)
adata = adata[adata.obs['n_genes_by_counts'] > 200, :]
adata = adata[adata.obs['pct_counts_mt'] < 10, :]

# === 5. Create pseudobulk counts ===
# Sum counts across all cells from same sample
pseudobulk_counts = (
    adata.to_df()
    .groupby(adata.obs['sample'])
    .sum()
    .T  # Transpose → genes × samples
)

# Save pseudobulk count matrix
pseudobulk_counts.to_csv("pseudobulk_counts.csv")

# Save sample metadata
meta = pd.DataFrame({
    'sample': ['control_1', 'treated_1'],
    'condition': ['control', 'treated']
})
meta.to_csv("sample_metadata.csv", index=False)