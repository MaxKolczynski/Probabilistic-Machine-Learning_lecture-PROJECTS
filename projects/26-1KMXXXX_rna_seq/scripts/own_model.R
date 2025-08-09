# 1. Libraries
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(MASS)

library(pbapply)
install.packages("pbapply")

# 2. Arbeitsverzeichnis
setwd("/Users/maxvk/Probabilistic-Machine-Learning_lecture-PROJECTS/projects/26-1KMXXXX_rna_seq")

# 3. Rohdaten einlesen
counts <- read.delim("data/GSE126848_raw_counts_GRCh38.p13_NCBI.tsv", 
                     header = TRUE, 
                     row.names = 1, 
                     check.names = FALSE)

# 4. Sample-IDs definieren
nash_samples <- c("GSM3615308", "GSM3615309", "GSM3615310", "GSM3615311", "GSM3615312",
                  "GSM3615313", "GSM3615314", "GSM3615315", "GSM3615316", "GSM3615317",
                  "GSM3615318", "GSM3615319", "GSM3615320", "GSM3615321", "GSM3615322")

healthy_samples <- c("GSM3615323", "GSM3615324", "GSM3615325", "GSM3615326", "GSM3615327",
                     "GSM3615328", "GSM3615329", "GSM3615330", "GSM3615331", "GSM3615332",
                     "GSM3615333", "GSM3615334", "GSM3615335", "GSM3615336")

selected_samples <- c(nash_samples, healthy_samples)

# 5. Subset der Zählmatrix auf ausgewählte Samples
counts_subset <- counts[, selected_samples]

# 6. Sample-Metadaten
condition <- factor(c(rep("NASH", length(nash_samples)), rep("Healthy", length(healthy_samples))),
                    levels = c("Healthy", "NASH"))

colData <- data.frame(
  sample = selected_samples,
  condition = condition
)

# 7. Auswahl eines real existierenden Gens mit hoher Gesamtanzahl
gene_sums <- rowSums(counts_subset)
gene_id <- names(sort(gene_sums, decreasing = TRUE))[1]  # Top-expressed Gene

# 8. Extrahiere Counts für das gewählte Gen
gene_counts <- counts_subset[gene_id, ]

# 9. Kombiniere Counts mit Metadaten
df <- data.frame(
  counts = as.numeric(gene_counts),
  sample = names(gene_counts)
) %>%
  left_join(colData, by = "sample")

# 10. Sicherstellen, dass Condition ein Faktor ist
df$condition <- factor(df$condition)

# 11. NB-GLM fitten
fit <- glm.nb(counts ~ condition, data = df)
summary(fit)


# Liste aller Gene
gene_list <- rownames(counts_subset)

# Ergebnisse speichern
glm_results <- lapply(gene_list, function(gene_id) {
  gene_counts <- counts_subset[gene_id, ]
  
  df <- data.frame(
    counts = as.numeric(gene_counts),
    sample = names(gene_counts)
  ) %>%
    left_join(colData, by = "sample")
  
  df$condition <- factor(df$condition)
  
  # Modell fitten mit Fehlerbehandlung
  tryCatch({
    fit <- glm.nb(counts ~ condition, data = df)
    coefs <- summary(fit)$coefficients
    
    data.frame(
      gene = gene_id,
      logFC = coefs["conditionNASH", "Estimate"],
      se = coefs["conditionNASH", "Std. Error"],
      pval = coefs["conditionNASH", "Pr(>|z|)"]
    )
  }, error = function(e) {
    # Wenn das Modell nicht konvergiert
    data.frame(
      gene = gene_id,
      logFC = NA,
      se = NA,
      pval = NA
    )
  })
})

# In DataFrame umwandeln
results_df <- bind_rows(glm_results)

# Multiple Testing Correction
results_df$padj <- p.adjust(results_df$pval, method = "BH")




# Volcano Plot
ggplot(results_df, aes(x = logFC, y = -log10(padj))) +
  geom_point(alpha = 0.4) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  labs(title = "Volcano Plot", x = "log2 Fold Change", y = "-log10 Adjusted p-value") +
  theme_minimal()

results_df <- results_df %>%
  mutate(significant = ifelse(padj < 0.05 & abs(logFC) > 1, "yes", "no"))

ggplot(results_df, aes(x = logFC, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("no" = "grey", "yes" = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "Volcano Plot mit Signifikanz", x = "log2 Fold Change", y = "-log10 Adjusted p-value") +
  theme_minimal()


# other version 
library(ggrepel)

# Setze log2FC-Schwelle und p-Adjust-Schwelle
log2fc_cutoff <- 1
padj_cutoff <- 0.05

# Bereinige DataFrame
res_df <- results_df %>%
  filter(!is.na(padj) & !is.na(logFC)) %>%
  rename(log2FoldChange = logFC)

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


# other version - final 
library(ggrepel)

# Setze log2FC-Schwelle und p-Adjust-Schwelle
log2fc_cutoff <- 1
padj_cutoff <- 0.05

# Bereinige DataFrame
res_df <- results_df %>%
  filter(!is.na(padj) & !is.na(logFC)) %>%
  rename(log2FoldChange = logFC)

# Klassifiziere Gene
res_df$significance <- "ns"
res_df$significance[res_df$log2FoldChange >  log2fc_cutoff & res_df$padj < padj_cutoff] <- "up"
res_df$significance[res_df$log2FoldChange < -log2fc_cutoff & res_df$padj < padj_cutoff] <- "down"

# Optional: Top Gene labeln
top_genes <- res_df %>%
  filter(padj < 1e-10 & abs(log2FoldChange) > 2)

# Volcano Plot mit festen Achsgrenzen
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_text_repel(data = top_genes, aes(label = gene), size = 3, max.overlaps = 20) +
  scale_color_manual(values = c("up" = "red", "down" = "blue", "ns" = "grey")) +
  geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed") +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
  coord_cartesian(xlim = c(-8, 8), ylim = c(0, 35)) +
  theme_minimal() +
  labs(title = "Volcano Plot: NASH vs Healthy",
       x = "log2 Fold Change",
       y = "-log10 Adjusted p-value",
       color = "Regulation")



# other version
library(ggrepel)

# Ergebnisse als DataFrame
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# Entferne NA-Werte
res_df <- na.omit(res_df)

# Signifikanzschwellen
log2fc_cutoff <- 1
padj_cutoff <- 0.05

# Bereinige DataFrame
res_df <- results_df %>%
  filter(!is.na(padj) & !is.na(logFC)) %>%
  rename(log2FoldChange = logFC)

# Klassifiziere Gene
res_df$significance <- "ns"
res_df$significance[res_df$log2FoldChange >  log2fc_cutoff & res_df$padj < padj_cutoff] <- "up"
res_df$significance[res_df$log2FoldChange < -log2fc_cutoff & res_df$padj < padj_cutoff] <- "down"

# Optional: Top-Gene labeln
top_genes <- res_df %>%
  filter(padj < 1e-10 & abs(log2FoldChange) > 2)

# Dynamische Achsengrenzen berechnen
xlim_max <- ceiling(max(abs(res_df$log2FoldChange), na.rm = TRUE)) + 0.5
ylim_max <- ceiling(max(-log10(res_df$padj), na.rm = TRUE)) + 1

# Volcano Plot
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_text_repel(data = top_genes, aes(label = gene), size = 3, max.overlaps = 20) +
  scale_color_manual(values = c("up" = "red", "down" = "blue", "ns" = "grey")) +
  geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed") +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
  coord_cartesian(xlim = c(-xlim_max, xlim_max), ylim = c(0, ylim_max)) +
  theme_minimal() +
  labs(title = "Volcano Plot: NASH vs Healthy",
       x = "log2 Fold Change",
       y = "-log10 Adjusted p-value",
       color = "Regulation")




write.csv(results_df, file = "nb_glm_differential_expression_results.csv", row.names = FALSE)

############ END

# Schritt 1: Pseudocounts vermeiden
counts_pseudo <- counts_subset + 1

# Schritt 2: Berechne geometrisches Mittel je Gen
geo_means <- apply(counts_pseudo, 1, function(x) exp(mean(log(x))))

# Schritt 3: Ratio jedes Counts zur geom. Mitte
ratios <- sweep(counts_pseudo, 1, geo_means, FUN = "/")

# Schritt 4: Median Ratio pro Sample → Size Factors
size_factors <- apply(ratios, 2, median, na.rm = TRUE)

# Normalize Counts (wie counts(dds, normalized=TRUE))
norm_counts <- sweep(counts_subset, 2, size_factors, FUN = "/")

# Berechne Mittelwert und Varianz pro Gen
mu <- rowMeans(norm_counts)
var <- apply(norm_counts, 1, var)

# Roh-Dispersion pro Gen nach Formel: Var(Y) = mu + alpha * mu^2 → alpha = (Var - mu) / mu^2
disp_raw <- (var - mu) / mu^2
disp_raw[disp_raw < 0] <- NA  # Filter unplausible Werte

disp_df <- data.frame(
  mean = mu,
  dispersion = disp_raw
)

ggplot(disp_df, aes(x = mean, y = dispersion)) +
  geom_point(alpha = 0.4) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = "Dispersion Estimates (raw)", x = "Mean normalized count", y = "Dispersion") +
  theme_minimal()

dds <- DESeqDataSetFromMatrix(countData = counts_subset,
                              colData = colData,
                              design = ~ condition)

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
plotDispEsts(dds)


deg_table <- results_df %>%
  filter(!is.na(logFC) & !is.na(padj)) %>%
  rename(
    log2FoldChange = logFC,
    std_error = se
  ) %>%
  dplyr::select(gene, log2FoldChange, std_error, pval, padj) %>%
  arrange(padj)

head(deg_table, 10)



