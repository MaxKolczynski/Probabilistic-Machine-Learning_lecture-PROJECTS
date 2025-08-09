# ------------------ 1. 📚 Libraries ------------------
library(tidyverse)
library(MASS)
library(ggplot2)
library(ggrepel)
library(DESeq2)

# ------------------ 2. 📁 Arbeitsverzeichnis ------------------
setwd("/Users/maxvk/Probabilistic-Machine-Learning_lecture-PROJECTS/projects/26-1KMXXXX_rna_seq")

# ------------------ 3. 📄 Daten einlesen ------------------
counts <- read.delim("data/GSE126848_raw_counts_GRCh38.p13_NCBI.tsv",
                     header = TRUE,
                     row.names = 1,
                     check.names = FALSE)

# ------------------ 4. 🎯 Sample-Subset & Umbenennung via Metadaten ------------------

# A. Metadaten einlesen
metadata <- read.csv("data/metadata.csv", stringsAsFactors = FALSE, check.names = FALSE)

# B. Relevante Spalten umbenennen (achte auf Backticks!)
metadata <- metadata %>%
  rename(
    sample_id = `GEO_Accession (exp)`,
    condition_raw = disease
  )

# C. Nur Samples mit gewünschter Erkrankung behalten
metadata <- metadata %>%
  mutate(condition = case_when(
    grepl("NAFLD", condition_raw, ignore.case = TRUE) ~ "NASH",
    grepl("healthy|normal|control", condition_raw, ignore.case = TRUE) ~ "Healthy",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(condition), sample_id %in% colnames(counts))

# D. Neue Sample-Namen erzeugen (z. B. NASH_1, Healthy_2, …)
metadata <- metadata %>%
  group_by(condition) %>%
  mutate(new_name = paste0(condition, "_", row_number())) %>%
  ungroup()

# E. Count-Matrix subsetten und Spalten umbenennen
counts_subset <- counts[, metadata$sample_id]
colnames(counts_subset) <- metadata$new_name

# F. colData vorbereiten
colData <- metadata %>%
  select(sample = new_name, condition)
rownames(colData) <- colData$sample
colData$condition <- factor(colData$condition, levels = c("Healthy", "NASH"))

# A. Originale Sample-IDs
#nash_samples <- c("GSM3615308", "GSM3615309", "GSM3615310", "GSM3615311", "GSM3615312",
#                  "GSM3615313", "GSM3615314", "GSM3615315", "GSM3615316", "GSM3615317",
#                  "GSM3615318", "GSM3615319", "GSM3615320", "GSM3615321", "GSM3615322")

#healthy_samples <- c("GSM3615323", "GSM3615324", "GSM3615325", "GSM3615326", "GSM3615327",
#                     "GSM3615328", "GSM3615329", "GSM3615330", "GSM3615331", "GSM3615332",
#                     "GSM3615333", "GSM3615334", "GSM3615335", "GSM3615336")

# B. Neue Labels definieren
#nash_labels <- paste0("NASH_", seq_along(nash_samples))
#healthy_labels <- paste0("Healthy_", seq_along(healthy_samples))

# C. Sample-Mapping erstellen (alt → neu)
#sample_rename_map <- c(setNames(nash_labels, nash_samples),
                       setNames(healthy_labels, healthy_samples))

# D. Subset der Count-Matrix und Umbenennung der Spalten
#selected_samples <- c(nash_samples, healthy_samples)
#counts_subset <- counts[, selected_samples]
#colnames(counts_subset) <- sample_rename_map[colnames(counts_subset)]

# E. Metadaten-Frame (colData) mit neuen Namen
#condition <- factor(c(rep("NASH", length(nash_samples)),
#                      rep("Healthy", length(healthy_samples))),
#                    levels = c("Healthy", "NASH"))

#colData <- data.frame(
#  sample = sample_rename_map[selected_samples],
#  condition = condition
#)
#rownames(colData) <- colData$sample  # wichtig für DESeq2-Kompatibilität



# ------------------ 4b. QC nach HBC Training ------------------
# ------------------ 4b. QC: Sample-Level & Gene-Level QC (nach HBC-Tutorial) ------------------

# A. 📏 Library Size Check (Sequenziertiefe je Sample)
libsize <- colSums(counts_subset)
barplot(libsize, las = 2, col = "steelblue",
        main = "Library Sizes", ylab = "Total Counts")

# B. 🔍 Vergleich: klassische Mindestfilterung (zur Referenz)
keep_simple <- rowSums(counts_subset >= 10) >= 2
cat("Gene mit ≥10 Counts in ≥2 Samples (klassischer Filter):", sum(keep_simple), "\n")

# C. 🧬 Gene-Level QC – mehrstufig nach HBC
# 1. Entferne Gene mit nur 0 Counts in allen Proben
counts_nonzero <- counts_subset[rowSums(counts_subset) > 0, ]
cat("Gene übrig nach Entfernen von all-zero:", nrow(counts_nonzero), "\n")

# 2. Entferne Gene mit extremen Count-Ausreißern (MAD-basierter Filter)
mad_gene <- apply(counts_nonzero, 1, mad)
keep_mad <- mad_gene < quantile(mad_gene, 0.99)  # oberstes 1% wird entfernt
counts_mad_filtered <- counts_nonzero[keep_mad, ]
cat("Gene übrig nach Entfernen extremer Ausreißer:", nrow(counts_mad_filtered), "\n")

# 3. Entferne Gene mit geringer mittlerer Expression
mean_counts <- rowMeans(counts_mad_filtered)
keep_mean <- mean_counts > 5  # Schwelle je nach Datensatz anpassbar
counts_filtered <- counts_mad_filtered[keep_mean, ]
cat("Gene übrig nach Entfernen von low-expression-Genen:", nrow(counts_filtered), "\n")

# D. ✨ VST-Transformation (nur für Visualisierung, nicht für Modell)
dds_vst <- DESeqDataSetFromMatrix(countData = counts_filtered,
                                  colData = colData,
                                  design = ~ condition)
dds_vst <- estimateSizeFactors(dds_vst)
vsd <- vst(dds_vst, blind = FALSE)

# E. 📊 PCA zur QC (Sample-Level)
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% Var")) +
  ylab(paste0("PC2: ", percentVar[2], "% Var")) +
  theme_minimal() +
  ggtitle("PCA: Sample-Level QC")

# F. 🌡️ Sample Correlation Heatmap + Hierarchical Clustering
logcounts <- assay(vsd)
sample_cor <- cor(logcounts)
pheatmap::pheatmap(sample_cor,
                   main = "Sample Correlation Heatmap",
                   display_numbers = TRUE,
                   color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(100),
                   clustering_method = "complete")

# G. 📦 Boxplot der log-normalisierten Counts
boxplot(logcounts, las = 2, col = "lightgray",
        main = "Log-Counts (VST)", ylab = "log2(count + 1)")


# ------------------ 5. Normalisierung (manuell) ------------------
counts_pseudo <- counts_subset + 1
geo_means <- apply(counts_pseudo, 1, function(x) exp(mean(log(x))))
ratios <- sweep(counts_pseudo, 1, geo_means, FUN = "/")
size_factors <- apply(ratios, 2, median, na.rm = TRUE)
norm_counts <- sweep(counts_subset, 2, size_factors, FUN = "/")

# ------------------ 6. Dispersion Plot (DESeq2-Style) ------------------

# Erstelle DESeqDataSet-Objekt
dds <- DESeqDataSetFromMatrix(countData = counts_subset,
                              colData = colData,
                              design = ~ condition)

# Schätze Size Factors und Dispersions (nur notwendig für Plot)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

# DESeq2-eigener Plot der genauen + gefitteten Dispersionsschätzungen
plotDispEsts(dds)

# ------------------ 7. GLM: Gene-wise NB-Fits (mit Progress Bar) ------------------

# Fortschrittsanzeige aktivieren
library(pbapply)

gene_list <- rownames(counts_subset)

# pblapply statt lapply → zeigt Fortschrittsbalken
glm_results <- pblapply(gene_list, function(gene_id) {
  gene_counts <- counts_subset[gene_id, ]
  df <- data.frame(counts = as.numeric(gene_counts),
                   sample = names(gene_counts)) %>%
    left_join(colData, by = "sample")
  df$condition <- factor(df$condition)
  
  tryCatch({
    fit <- glm.nb(counts ~ condition, data = df, link = log)
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

# Kombiniere Ergebnisse
results_df <- bind_rows(glm_results)
results_df$padj <- p.adjust(results_df$pval, method = "BH")

# ------------------ 8. Volcano Plot – FINAL VERSION ------------------
log2fc_cutoff <- 1
padj_cutoff <- 0.05

res_df <- results_df %>%
  filter(!is.na(padj) & !is.na(logFC)) %>%
  rename(log2FoldChange = logFC)

res_df$significance <- "ns"
res_df$significance[res_df$log2FoldChange >  log2fc_cutoff & res_df$padj < padj_cutoff] <- "up"
res_df$significance[res_df$log2FoldChange < -log2fc_cutoff & res_df$padj < padj_cutoff] <- "down"

top_genes <- res_df %>%
  filter(padj < 1e-10 & abs(log2FoldChange) > 2)

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_text_repel(data = top_genes, aes(label = gene), size = 3, max.overlaps = 20) +
  scale_color_manual(values = c("up" = "red", "down" = "blue", "ns" = "grey")) +
  geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed") +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
  coord_cartesian(xlim = c(-7, 7), ylim = c(0, 35)) +
  theme_minimal() +
  labs(title = "Volcano Plot: NASH vs Healthy",
       x = "log2 Fold Change", y = "-log10 Adjusted p-value", color = "Regulation")

# ------------------ 9. DEG-Tabelle exportieren ------------------
write.csv(results_df, file = "nb_glm_differential_expression_results.csv", row.names = FALSE)

# ------------------ 10. Top DEGs anzeigen ------------------
# Berechne baseMean über alle normalisierten Counts
base_means <- rowMeans(norm_counts)

# Erstelle DEG-Tabelle im DESeq2-Stil
deg_table <- results_df %>%
  filter(!is.na(logFC) & !is.na(pval) & !is.na(padj)) %>%
  mutate(
    baseMean = base_means[gene],
    stat = logFC / se  # z-Wert: Effekt / Standardfehler
  ) %>%
  rename(
    log2FoldChange = logFC,
    lfcSE = se,
    pvalue = pval
  ) %>%
  dplyr::select(gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>%
  arrange(padj)

# Zeige die Top 10
head(deg_table, 10)
