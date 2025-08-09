# ------------------ 1. 📚 Libraries ------------------
library(tidyverse)
library(MASS)
library(ggplot2)
library(ggrepel)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(pbapply)

# ------------------ 2. 📁 Arbeitsverzeichnis ------------------
setwd("/Users/maxvk/Probabilistic-Machine-Learning_lecture-PROJECTS/projects/26-1KMXXXX_rna_seq")

# ------------------ 3. 📄 Daten einlesen ------------------
counts <- read.delim("data/GSE126848_raw_counts_GRCh38.p13_NCBI.tsv",
                     header = TRUE,
                     row.names = 1,
                     check.names = FALSE)

metadata <- read.csv("data/metadata.csv", stringsAsFactors = FALSE, check.names = FALSE)

# ------------------ 4. 🧪 Sample-Subset & Umbenennung ------------------

metadata <- metadata %>%
  rename(
    sample_id = `GEO_Accession (exp)`,
    condition_raw = disease
  ) %>%
  mutate(condition = case_when(
    grepl("NAFLD", condition_raw, ignore.case = TRUE) ~ "NASH",
    grepl("healthy|normal|control", condition_raw, ignore.case = TRUE) ~ "Healthy",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(condition), sample_id %in% colnames(counts)) %>%
  group_by(condition) %>%
  mutate(new_name = paste0(condition, "_", row_number())) %>%
  ungroup()

counts_subset <- counts[, metadata$sample_id]
colnames(counts_subset) <- metadata$new_name


colData <- metadata %>%
  dplyr::select(sample = new_name, condition)

rownames(colData) <- colData$sample
colData$condition <- factor(colData$condition, levels = c("Healthy", "NASH"))

# dataset summary
nrow(counts_subset)  # Anzahl der Gene (Zeilen)
summary(counts_subset)  # Zeigt Min, Median, Mean etc. über alle Werte
library(matrixStats)
rowSumsZero <- sum(rowSums(counts_subset) == 0)  # Gene mit 0 Reads über alle Samples
rowMeans <- rowMeans(counts_subset)              # Durchschnittliche Expression pro Gen
hist(log10(rowMeans + 1), main = "Log10 Mean Expression per Gene")  # Verteilung


# ------------------ 5. 🔍 Sample- & Gene-Level QC ------------------

# 📏 Library Sizes
barplot(colSums(counts_subset), las = 2, col = "steelblue", main = "Library Sizes", ylab = "Total Counts")

# Gene QC
keep_zero <- rowSums(counts_subset) > 0
counts_nonzero <- counts_subset[keep_zero, ]

mad_vals <- apply(counts_nonzero, 1, mad)
keep_mad <- mad_vals < quantile(mad_vals, 0.99)
counts_mad_filtered <- counts_nonzero[keep_mad, ]

mean_counts <- rowMeans(counts_mad_filtered)
keep_mean <- mean_counts > 5
counts_filtered <- round(as.matrix(counts_mad_filtered[keep_mean, ]))  # <-- wichtig

# ------------------ 6. ✨ VST + Visual QC ------------------

dds_vst <- DESeqDataSetFromMatrix(countData = counts_filtered,
                                  colData = colData,
                                  design = ~ condition)
dds_vst <- estimateSizeFactors(dds_vst)
vsd <- vst(dds_vst, blind = FALSE)



# PCA
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  ggtitle("PCA: Sample-Level QC") +
  theme_minimal()


# 1. Log2(count + 1)-Transformation auf normalisierte Counts
logcounts <- log2(counts(dds, normalized = TRUE) + 1)

# 2. PCA berechnen (Samples müssen in Zeilen stehen → Transponieren)
pca <- prcomp(t(logcounts))

# 3. PCA-Ergebnisse mit Sample-Metadaten verbinden
pcaData <- as.data.frame(pca$x)
pcaData$condition <- colData(dds)$condition  # ggf. anpassen, z. B. "batch", "group", etc.

# 4. Prozentuale Varianz je PC berechnen
percentVar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)))

# 5. PCA-Plot mit ggplot2
library(ggplot2)
ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 6) +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  ggtitle("PCA: log2(count + 1) transformation") +
  theme_minimal()


# Heatmap
# Hole sortierte Sample-Namen nach Bedingung
sample_order <- colData %>%
  arrange(condition) %>%
  pull(sample)

# Sortiere die Korrelationsmatrix entsprechend
sample_cor_ordered <- sample_cor[sample_order, sample_order]

# Zeichne Heatmap ohne Clustering (bereits sortiert)
pheatmap::pheatmap(sample_cor_ordered,
                   main = "Sample Correlation Heatmap (ordered by condition)",
                   display_numbers = TRUE,
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(100))
# Boxplot
boxplot(logcounts, las = 2, col = "lightgray",
        main = "Log-Counts", ylab = "log2(count + 1)")



# ------------------ 7. 📈 Dispersion Plot ------------------

dds <- DESeqDataSetFromMatrix(countData = counts_subset,
                              colData = colData,
                              design = ~ condition)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
plotDispEsts(dds)


# alternative 
library(DESeq2)
library(ggplot2)

# Set up DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_subset,
                              colData = colData,
                              design = ~ condition)

# Only estimate size factors and dispersion
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

# Extract mean normalized counts
mean_counts <- rowMeans(counts(dds, normalized = TRUE))

# Extract dispersion estimates
dispersions <- dispersions(dds)

# Create data frame
disp_df <- data.frame(
  mean = mean_counts,
  dispersion = dispersions
)

# Plot: Dispersion vs. Mean (log-log scale)
ggplot(disp_df, aes(x = mean, y = dispersion)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Gene-wise Dispersion Estimates (no shrinkage)",
    x = "Mean of normalized counts (log scale)",
    y = "Dispersion estimate (log scale)"
  ) +
  theme_minimal()

# alternative 
library(DESeq2)
library(ggplot2)

# Daten vorbereiten
dds <- DESeqDataSetFromMatrix(countData = counts_subset,
                              colData = colData,
                              design = ~ condition)
dds <- estimateSizeFactors(dds)

# Normalisierte Counts extrahieren
norm_counts <- counts(dds, normalized = TRUE)

# Mittelwert und Varianz berechnen (pro Gen = Zeile)
mean_counts <- rowMeans(norm_counts)
var_counts <- apply(norm_counts, 1, var)

# Dataframe für Plot
var_df <- data.frame(
  mean = mean_counts,
  variance = var_counts
)

# Plot erstellen
ggplot(var_df, aes(x = mean, y = variance)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Empirical Variance vs. Mean (Normalized Counts)",
    subtitle = "Red dashed line shows Poisson assumption (Var = Mean)",
    x = "Mean of normalized counts (log10)",
    y = "Empirical variance (log10)"
  ) +
  theme_minimal()

# ------------------ 8. 🔧 Manuelle Normalisierung für GLM ------------------

# 1. Berechnung der size factors (wie gehabt)
counts_pseudo <- counts_subset + 1
geo_means <- apply(counts_pseudo, 1, function(x) exp(mean(log(x))))
ratios <- sweep(counts_pseudo, 1, geo_means, FUN = "/")
size_factors <- apply(ratios, 2, median, na.rm = TRUE)

# 2. Sample-Zuordnung der size factors
sample_size_factors <- size_factors
names(sample_size_factors) <- colnames(counts_subset)

# ------------------ 9. 🔬 GLM: Negative Binomial für jedes Gen ------------------


#gene_list <- rownames(counts_subset)
#glm_results <- pblapply(gene_list, function(gene_id) {
#  gene_counts <- counts_subset[gene_id, ]
#  df <- data.frame(counts = as.numeric(gene_counts),
#                   sample = names(gene_counts)) %>%
#    left_join(colData, by = "sample")
#  df$condition <- factor(df$condition)
#  
#  tryCatch({
#    fit <- glm.nb(counts ~ condition, data = df, link = log)
#    coefs <- summary(fit)$coefficients
#    data.frame(
#      gene = gene_id,
#      logFC = coefs["conditionNASH", "Estimate"],
#      se = coefs["conditionNASH", "Std. Error"],
#      pval = coefs["conditionNASH", "Pr(>|z|)"]
#    )
#  }, error = function(e) {
#    data.frame(gene = gene_id, logFC = NA, se = NA, pval = NA)
#  })
#})

#results_df <- bind_rows(glm_results)
#results_df$padj <- p.adjust(results_df$pval, method = "BH")



# 3. GLM-Schleife über Gene mit Offset
gene_list <- rownames(counts_subset)
glm_results <- pblapply(gene_list, function(gene_id) {
  gene_counts <- counts_subset[gene_id, ]
  
  df <- data.frame(
    counts = as.numeric(gene_counts),
    sample = names(gene_counts)
  ) %>%
    left_join(colData, by = "sample")
  
  df$condition <- relevel(factor(df$condition), ref = "Healthy")
  df$size_factor <- sample_size_factors[df$sample]
  
  tryCatch({
    fit <- glm.nb(counts ~ condition + offset(log(size_factor)), data = df, link = log)
    coefs <- summary(fit)$coefficients
    data.frame(
      gene = gene_id,
      logFC = coefs["conditionNASH", "Estimate"],
      se = coefs["conditionNASH", "Std. Error"],
      pval = coefs["conditionNASH", "Pr(>|z|)"]
    )
  }, error = function(e) {
    data.frame(gene = gene_id, logFC = NA, se = NA, pval = NA)
  })
})

# 4. Ergebnis-Dataframe + multiple testing correction
results_df <- bind_rows(glm_results)
results_df$padj <- p.adjust(results_df$pval, method = "BH")

sum(results_df$padj < 0.05, na.rm = TRUE)

# Anzahl hochregulierter Gene (logFC > 0)
sum(results_df$padj < 0.05 & results_df$logFC > 0, na.rm = TRUE)

# Anzahl runterregulierter Gene (logFC < 0)
sum(results_df$padj < 0.05 & results_df$logFC < 0, na.rm = TRUE)

# ------------------ 10. 🌋 Volcano Plot ------------------

log2fc_cutoff <- 1
padj_cutoff <- 0.05

res_df <- results_df %>%
  filter(!is.na(padj) & !is.na(logFC)) %>%
  rename(log2FoldChange = logFC) %>%
  mutate(
    significance = case_when(
      log2FoldChange > log2fc_cutoff & padj < padj_cutoff ~ "up",
      log2FoldChange < -log2fc_cutoff & padj < padj_cutoff ~ "down",
      TRUE ~ "ns"
    )
  )

top_genes <- res_df %>%
  filter(padj < 1e-10 & abs(log2FoldChange) > 2)

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_text_repel(data = top_genes, aes(label = gene), size = 3, max.overlaps = 20) +
  scale_color_manual(values = c("up" = "red", "down" = "blue", "ns" = "grey")) +
  geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed") +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
  coord_cartesian(xlim = c(-7, 7), ylim = c(0, 50)) +
  theme_minimal() +
  labs(title = "Volcano Plot: NASH vs Healthy",
       x = "log2 Fold Change", y = "-log10 Adjusted p-value", color = "Regulation")

# ------------------ 11. 💾 DEG-Tabelle speichern ------------------

write.csv(results_df, file = "nb_glm_differential_expression_results.csv", row.names = FALSE)

summary(results_df)

# ------------------ 12. 🧬 Top DEGs anzeigen ------------------

base_means <- rowMeans(norm_counts)

deg_table <- results_df %>%
  filter(!is.na(logFC) & !is.na(pval) & !is.na(padj)) %>%
  mutate(
    baseMean = base_means[gene],
    stat = logFC / se
  ) %>%
  rename(
    log2FoldChange = logFC,
    lfcSE = se,
    pvalue = pval
  ) %>%
  select(gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>%
  arrange(padj)

head(deg_table, 10)

