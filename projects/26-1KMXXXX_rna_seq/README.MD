# NB-GLM bulk RNA-seq (NASH vs. Healthy)

## Struktur

project/
├─ data/
│  ├─ GSE126848_raw_counts_GRCh38.p13_NCBI.tsv
│  └─ sample_lists.yaml
├─ outputs/
│  ├─ figures/
│  └─ tables/
├─ R/
│  ├─ 00_utils.R
│  ├─ 01_load_data.R
│  ├─ 02_filter_qc.R
│  ├─ 03_normalize.R
│  ├─ 04_qc_plots.R
│  ├─ 05_nb_glm.R
│  └─ 06_export.R
└─ main.R

## Setup

```r
install.packages(c(
  "tidyverse","here","yaml","matrixStats","MASS","pbapply",
  "pheatmap","RColorBrewer","ggrepel","conflicted"
))

