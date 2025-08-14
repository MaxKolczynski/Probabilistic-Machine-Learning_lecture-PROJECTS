
# 

This repository contains a probabilistic machine learning project for detecting differentially expressed genes using generalized linear models. 

# 🧬 Dataset 

The dataset is publicy available under the Gene Expression Omnibus (GEO) Acession number: [GSE126848](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126848)

The original Pubilcation can be found [here](https://pubmed.ncbi.nlm.nih.gov/30653341/)

# 📂 Repository Tree Structure

```
project/
├─ data/  
│  ├─ GSE126848_raw_counts_GRCh38.p13_NCBI.tsv   # count matrix
│  └─ sample_lists.yaml                          # sample mapping 
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
```

# Notes on the code modules

## 00_utils.R  
## 01_load_data.R  
## 02_filter_qc.R  
## 03_normalize.R  
## 04_qc_plots.R  
## 05_nb_glm.R  
## 06_export.R  

# Setup

Simply execute ```Rscript main.R`` to get the entire piepline running. This will automatically procude figures and tables that can be found in the ```outputs``` fodler. 

```r
install.packages(c(
  "tidyverse","here","yaml","matrixStats","MASS","pbapply",
  "pheatmap","RColorBrewer","ggrepel","conflicted"
))



