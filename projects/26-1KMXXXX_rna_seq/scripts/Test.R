
# 1. Libraries
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(MASS)

setwd("/Users/maxvk/Probabilistic-Machine-Learning_lecture-PROJECTS/projects/26-1KMXXXX_rna_seq")

counts <- read.delim("data/GSE126848_raw_counts_GRCh38.p13_NCBI.tsv", 
                     header = TRUE, 
                     row.names = 1, 
                     check.names = FALSE)

dim(counts)                   # Zeilen = Gene, Spalten = Samples?
head(rownames(counts), 5)     # Sind das Ensembl-IDs oder Gen-Namen?
colnames(counts)[1:10]        # Welche Sample-Namen gibt es?
str(counts[, 1:3])            # Sind es numerische Count-Werte?

colnames(counts)

# Vektor der GSM-IDs
nash_samples <- c("GSM3615308", "GSM3615309", "GSM3615310", "GSM3615311", "GSM3615312",
                  "GSM3615313", "GSM3615314", "GSM3615315", "GSM3615316", "GSM3615317",
                  "GSM3615318", "GSM3615319", "GSM3615320", "GSM3615321", "GSM3615322")

healthy_samples <- c("GSM3615323", "GSM3615324", "GSM3615325", "GSM3615326", "GSM3615327",
                     "GSM3615328", "GSM3615329", "GSM3615330", "GSM3615331", "GSM3615332",
                     "GSM3615333", "GSM3615334", "GSM3615335", "GSM3615336")

selected_samples <- c(nash_samples, healthy_samples)

# Subset counts
counts_subset <- counts[, selected_samples]

# Faktor definieren
condition <- factor(c(rep("NASH", length(nash_samples)), rep("Healthy", length(healthy_samples))),
                    levels = c("Healthy", "NASH"))

colData <- data.frame(row.names = selected_samples,
                      condition = condition)

#--------------------- my model --------------------- 

gene_id <- "GENE_ID"  # Ersetze mit einem tatsächlichen Gen, z. B. "ENSG00000130203"
gene_counts <- counts_subset[gene_id, ]

df <- data.frame(
  counts = as.numeric(gene_counts),
  sample = names(gene_counts)
)

# Merge mit colData
df <- left_join(df, 
                tibble::rownames_to_column(colData, var = "sample"), 
                by = "sample")

# Check Levels
df$condition <- factor(df$condition)  # Erzwinge Faktor

# Modell fitten
fit <- glm.nb(counts ~ condition, data = df)
summary(fit)


results <- lapply(rownames(counts_subset), function(gene) {
  y <- as.numeric(counts_subset[gene, ])
  tryCatch({
    model <- glm.nb(y ~ condition, data = colData)
    coef_summary <- summary(model)$coefficients
    data.frame(
      gene = gene,
      logFC = coef_summary["conditionNASH", "Estimate"],
      se = coef_summary["conditionNASH", "Std. Error"],
      pval = coef_summary["conditionNASH", "Pr(>|z|)"]
    )
  }, error = function(e) {
    data.frame(
      gene = gene,
      logFC = NA,
      se = NA,
      pval = NA
    )
  })
})

results_df <- bind_rows(results)

results_df$padj <- p.adjust(results_df$pval, method = "BH")

ggplot(results_df, aes(x = logFC, y = -log10(padj))) +
  geom_point(alpha = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal()



#--------------------- DESeq2 --------------------- 
library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = counts_subset,
                              colData = colData,
                              design = ~ condition)

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
plotDispEsts(dds)

# Optionale Filterung: Gene mit sehr niedrigen Counts entfernen
dds <- dds[rowSums(counts(dds)) >= 10, ]


# Starte DESeq2
dds <- DESeq(dds)

# Ergebnisse extrahieren
res <- results(dds)

# Ergebnisse checken
summary(res)


library(ggplot2)
library(ggrepel)

# Ergebnisse als DataFrame
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# Entferne NA-Werte
res_df <- na.omit(res_df)

# Signifikanzschwellen
log2fc_cutoff <- 1
padj_cutoff <- 0.05

# Klassifiziere Gene
res_df$significance <- "ns"
res_df$significance[res_df$log2FoldChange >  log2fc_cutoff & res_df$padj < padj_cutoff] <- "up"
res_df$significance[res_df$log2FoldChange < -log2fc_cutoff & res_df$padj < padj_cutoff] <- "down"

# Optional: Top Gene labeln
top_genes <- res_df[res_df$padj < 1e-10 & abs(res_df$log2FoldChange) > 2, ]

# Volcano Plot
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_text_repel(data = top_genes, aes(label = gene), size = 3, max.overlaps = 20) +
  scale_color_manual(values = c("up" = "red", "down" = "blue", "ns" = "grey")) +
  geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed") +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
  theme_minimal() +
  labs(title = "Volcano Plot: NASH vs Healthy",
       x = "log2 Fold Change",
       y = "-log10 Adjusted p-value",
       color = "Regulation")


