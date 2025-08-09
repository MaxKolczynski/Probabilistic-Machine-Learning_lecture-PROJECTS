# 1. Libraries
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(pheatmap)

setwd("/Users/maxvk/Probabilistic-Machine-Learning_lecture-PROJECTS/projects/26-1KMXXXX_rna_seq") 

# Alle .txt.gz Dateien im RAW-Ordner listen
files <- list.files("data/GSE162045_RAW", pattern="*.txt.gz", full.names=TRUE)

# Dateinamen anzeigen
print(files)

# Schritt 2: Sample-Infos für DESeq2 anlegen

sample_names <- c("PC9_1", "PC9_2", "PC9_3", "PC9_4", "PC9_5", "PC9_6")
treatment <- c("untreated", "untreated", "untreated", "treated", "treated", "treated")

sampleTable <- data.frame(
  sampleName = sample_names,
  fileName = files,
  treatment = treatment,
  stringsAsFactors = FALSE
)
rownames(sampleTable) <- sampleTable$sampleName

# Sample-Tabelle anzeigen
print(sampleTable)


library(readr)

# Zeige die ersten 10 Zeilen der ersten Datei
read_tsv(sampleTable$fileName[1], n_max=10)


