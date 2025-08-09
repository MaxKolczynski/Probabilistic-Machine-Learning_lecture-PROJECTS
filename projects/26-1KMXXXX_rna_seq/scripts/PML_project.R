
library(Seurat)
library(dplyr)
library(ggplot2)

# === 1. Data Loading ===
control_data <- Read10X(data.dir = 'data/GSE273941_RAW/GSM8440146') # Control
treated_data <- Read10X(data.dir = 'data/GSE273941_RAW/GSM8440147') # Treated

seurat_ctrl <- CreateSeuratObject(counts = control_data)

# Schau dir die Zellnamen an
head(colnames(seurat_ctrl))  # oder View(colnames(seurat_ctrl))

# Schau dir die Metadaten an
head(seurat_ctrl@meta.data)



