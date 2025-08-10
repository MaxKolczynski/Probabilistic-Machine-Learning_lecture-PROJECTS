# ================== RNA-seq GLM (Negative Binomial) — Clean, Linear, Commented ==================
# Goal: Data exploration + preprocessing → manual normalization → gene-wise NB-GLM with offset →
#       multiple testing correction → export results → top DEGs → volcano plot
# Note: DESeq2 is NOT used for modeling; only base R / tidyverse / MASS. Comments in English.

# ------------------ 1) Libraries ------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(MASS)         # glm.nb
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(pbapply)      # progress apply
  library(matrixStats)  # efficient row/col ops
})

# ------------------ 2) Working directory ------------------
# Adjust if needed
setwd("/Users/maxvk/Probabilistic-Machine-Learning_lecture-PROJECTS/projects/26-1KMXXXX_rna_seq")

# ------------------ 3) Load raw data ------------------
counts_raw <- read.delim("data/GSE126848_raw_counts_GRCh38.p13_NCBI.tsv",
                         header = TRUE, row.names = 1, check.names = FALSE)

cat("=== Dataset (raw) ===\n")
cat("Genes (rows):", nrow(counts_raw), "\n")
cat("Samples (cols):", ncol(counts_raw), "\n\n")

# --- Section 4: build colData (explicit IDs) ---
# Achtung: verwende überall *nash_samples* (lowercase), nicht NASH_samples
nash_samples <- c(
  "GSM3615308","GSM3615309","GSM3615310","GSM3615311","GSM3615312",
  "GSM3615313","GSM3615314","GSM3615315","GSM3615316","GSM3615317",
  "GSM3615318","GSM3615319","GSM3615320","GSM3615321","GSM3615322"
)
healthy_samples <- c(
  "GSM3615323","GSM3615324","GSM3615325","GSM3615326","GSM3615327",
  "GSM3615328","GSM3615329","GSM3615330","GSM3615331","GSM3615332",
  "GSM3615333","GSM3615334","GSM3615335","GSM3615336"
)
selected_samples <- c(nash_samples, healthy_samples)

missing_in_counts <- setdiff(selected_samples, colnames(counts_raw))
if (length(missing_in_counts) > 0) {
  warning("Samples missing in count matrix and will be dropped: ",
          paste(missing_in_counts, collapse = ", "))
}
present_samples <- intersect(selected_samples, colnames(counts_raw))
stopifnot(length(present_samples) > 0)

counts <- counts_raw[, present_samples, drop = FALSE]

colData <- tibble::tibble(
  sample_id = present_samples,
  condition = ifelse(present_samples %in% nash_samples, "NASH", "Healthy")
) %>%
  dplyr::arrange(condition, sample_id) %>%
  dplyr::group_by(condition) %>%
  dplyr::mutate(sample = paste0(condition, "_", dplyr::row_number())) %>%
  dplyr::ungroup() %>%
  dplyr::select(sample, condition, sample_id)

# counts-Spalten in die neue Reihenfolge bringen und umbenennen
counts <- counts[, colData$sample_id, drop = FALSE]
colnames(counts) <- colData$sample
rownames(colData) <- colData$sample
colData$condition <- factor(colData$condition, levels = c("Healthy", "NASH"))

# Report condition counts
cat("=== Sample counts by condition (kept) ===\n")
print(colData %>% count(condition) %>% as.data.frame())
cat("Total kept samples:", nrow(colData), "\n\n")

# ------------------ 5) Quick data exploration (pre-filter) ------------------
# Library sizes
lib_sizes <- colSums(counts)
cat("=== Library sizes (raw counts) ===\n")
print(summary(lib_sizes)); cat("\n")

# Genes with all-zero counts
n_zero_genes <- sum(rowSums(counts) == 0)
cat("Genes with zero counts across all samples:", n_zero_genes, "\n\n")

# Histogram of mean expression (log10 scale)
row_means_raw <- rowMeans(counts)
hist(log10(row_means_raw + 1),
     main = "Log10 Mean Expression per Gene (raw)",
     xlab = "log10(mean+1)")

# Library sizes barplot
barplot(lib_sizes, las = 2, col = "steelblue",
        main = "Library Sizes (raw counts)", ylab = "Total counts")

# ------------------ 6) Gene-level filtering (QC) ------------------
# Rationale:
# - Remove genes with zero total counts (uninformative)
# - Remove extreme-MAD outliers (top 1%) to stabilize visuals
# - Keep genes with sufficient mean counts (> 5) to avoid quasi-separation
keep_nonzero <- rowSums(counts) > 0
counts_nz <- counts[keep_nonzero, , drop = FALSE]

mad_vals <- matrixStats::rowMads(as.matrix(counts_nz))
keep_mad <- mad_vals < quantile(mad_vals, 0.99, na.rm = TRUE)
counts_mad <- counts_nz[keep_mad, , drop = FALSE]

keep_mean <- rowMeans(counts_mad) > 5
counts_filt <- round(as.matrix(counts_mad[keep_mean, , drop = FALSE]))

# Report filtering impact
cat("=== Filtering summary ===\n")
cat("Genes before:", nrow(counts), "\n")
cat(" - after nonzero:", nrow(counts_nz), "\n")
cat(" - after MAD filter (99th pct removed):", nrow(counts_mad), "\n")
cat(" - after mean > 5:", nrow(counts_filt), "\n\n")

# ------------------ 7) Manual normalization (Median-of-Ratios, manual DESeq-style) ------------------
# Steps:
# 1) Add pseudocount to avoid log(0)
# 2) Per-gene geometric mean across samples
# 3) Ratios = counts / geo_mean per gene
# 4) Size factor = sample-wise median of ratios
# 5) Normalized counts = counts / size_factor
counts_pseudo <- counts_filt + 1
geo_means <- apply(counts_pseudo, 1, function(x) exp(mean(log(x))))
ratios <- sweep(counts_pseudo, 1, geo_means, FUN = "/")
size_factors <- apply(ratios, 2, median, na.rm = TRUE)
names(size_factors) <- colnames(counts_filt)

# Keep for modeling offset
sample_size_factors <- size_factors

# Normalized counts (no rounding)
norm_counts <- sweep(counts_filt, 2, size_factors, "/")

# Log2 normalized counts for visualization
logcounts <- log2(norm_counts + 1)

# Report size factor summary
cat("=== Size factors (median-of-ratios, manual) ===\n")
print(round(size_factors, 3))
cat("Summary:\n"); print(summary(size_factors)); cat("\n")

# --- ALIGNMENT FIX: Design exakt an logcounts ausrichten ---
sn <- colnames(logcounts)
design_df <- as.data.frame(colData)

# Set-Check mit klarer Diagnose bei Abweichungen
if (!setequal(sn, design_df$sample)) {
  cat("Only in logcounts:\n"); print(setdiff(sn, design_df$sample))
  cat("Only in colData:\n");   print(setdiff(design_df$sample, sn))
  stop("Sample name mismatch between logcounts and colData. Fix lists above.")
}

# Reihenfolge angleichen
design_aligned <- design_df[match(sn, design_df$sample), , drop = FALSE]
stopifnot(all(design_aligned$sample == sn))

# Überschreibe colData mit der ausgerichteten Version
colData <- design_aligned

# --- Section 8: PCA (ohne left_join, farbsicher) ---
pca <- prcomp(t(logcounts), center = TRUE, scale. = FALSE)
percentVar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 2)

pca_df <- as.data.frame(pca$x)
pca_df$sample <- colnames(logcounts)                 # identisch zu sn
pca_df$condition <- factor(colData$condition, levels = c("Healthy", "NASH"))

ggplot(pca_df, aes(PC1, PC2, color = condition, label = sample)) +
  geom_point(size = 3) +
  labs(
    title = "PCA: log2(normalized counts + 1)",
    x = paste0("PC1: ", percentVar[1], "%"),
    y = paste0("PC2: ", percentVar[2], "%")
  ) +
  theme_minimal()

# Boxplot of logcounts
boxplot(logcounts, las = 2, col = "lightgray",
        main = "Log2 Normalized Counts", ylab = "log2(norm+1)")

# Sample correlation heatmap (ordered by condition)
sample_cor <- cor(logcounts, method = "pearson")
sample_order <- colData %>% dplyr::arrange(condition) %>% dplyr::pull(sample)
sample_cor_ordered <- sample_cor[sample_order, sample_order]
pheatmap::pheatmap(sample_cor_ordered,
                   main = "Sample Correlation Heatmap (ordered by condition)",
                   display_numbers = TRUE,
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(100))

# ------------------ 9) Overdispersion diagnostics (no DESeq2) ------------------
# Empirical mean-variance relationship on normalized counts; Poisson baseline: Var ≈ Mean
mean_nc <- rowMeans(norm_counts)
var_nc  <- matrixStats::rowVars(as.matrix(norm_counts))

var_df <- data.frame(mean = mean_nc, variance = var_nc) %>%
  filter(mean > 0)

#ggplot(var_df, aes(x = mean, y = variance)) +
#  geom_point(alpha = 0.4, size = 1) +
#  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
#  scale_x_log10() + scale_y_log10() +
#  labs(title = "Empirical Variance vs Mean (Normalized Counts)",
#       subtitle = "Dashed line: Poisson assumption (Var = Mean)",
#       x = "Mean (log10)", y = "Variance (log10)") +
#  theme_minimal()

ggplot(var_df, aes(x = mean, y = variance)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  scale_x_log10() + scale_y_log10() +
  labs(
    title = "Empirical Variance vs Mean (Normalized Counts)",
    subtitle = "Dashed line: Poisson assumption (Var = Mean)",
    x = "Mean (log10)",
    y = "Variance (log10)"
  ) +
  theme_minimal()

# Optional: crude NB dispersion (alpha) estimate per gene via method-of-moments:
# Var = mu + alpha * mu^2  →  alpha ≈ max( (Var - mu) / mu^2 , 0 )
alpha_mom <- pmax((var_nc - mean_nc) / (mean_nc^2 + 1e-12), 0)
alpha_df <- data.frame(mean = mean_nc, alpha = alpha_mom) %>% filter(mean > 0)

ggplot(alpha_df, aes(x = mean, y = alpha)) +
  geom_point(alpha = 0.4, size = 1) +
  scale_x_log10() + scale_y_log10() +
  labs(title = "Method-of-Moments NB Dispersion (alpha) vs Mean",
       x = "Mean of normalized counts (log10)",
       y = "NB dispersion alpha (log10)") +
  theme_minimal()

# Report dispersion summaries
cat("=== Overdispersion (MoM alpha) summary ===\n")
print(summary(alpha_mom[is.finite(alpha_mom)]))
cat("\n")

# ------------------ 10) Gene-wise Negative Binomial GLM with offset ------------------
# Model: counts_g ~ condition + offset(log(size_factor))
# - Reference level: Healthy
# - Coefficient of interest: conditionNASH (log fold-change on NB scale)
gene_list <- rownames(counts_filt)
pboptions(type = "txt", style = 3)

# Prepare a small design tibble for joins

design_tbl <- as.data.frame(colData)[, c("sample","condition")]

glm_results <- pblapply(gene_list, function(gene_id) {
  y <- counts_filt[gene_id, ]
  df <- tibble::tibble(
    counts = as.numeric(y),
    sample = names(y)
  ) %>%
    dplyr::left_join(design_tbl, by = "sample") %>%
    dplyr::mutate(
      condition = stats::relevel(factor(condition), ref = "Healthy"),
      size_factor = sample_size_factors[sample]
    )
  
  # Defensive checks
  if (any(is.na(df$condition)) || any(is.na(df$size_factor))) {
    return(tibble(gene = gene_id, logFC = NA_real_, se = NA_real_, pval = NA_real_))
  }
  
  # Fit NB-GLM; use offset for library-size normalization
  fit <- try(glm.nb(counts ~ condition + offset(log(size_factor)), data = df, link = log),
             silent = TRUE)
  if (inherits(fit, "try-error")) {
    return(tibble(gene = gene_id, logFC = NA_real_, se = NA_real_, pval = NA_real_))
  }
  
  coefs <- summary(fit)$coefficients
  # Some genes may drop the level due to separation; guard access
  if (!("conditionNASH" %in% rownames(coefs))) {
    return(tibble(gene = gene_id, logFC = NA_real_, se = NA_real_, pval = NA_real_))
  }
  
  tibble(
    gene = gene_id,
    logFC = coefs["conditionNASH", "Estimate"],
    se    = coefs["conditionNASH", "Std. Error"],
    pval  = coefs["conditionNASH", "Pr(>|z|)"]
  )
})



results_df <- bind_rows(glm_results) %>%
  mutate(padj = p.adjust(pval, method = "BH"))

# Quick DEG counts (for report)
n_sig <- sum(results_df$padj < 0.05, na.rm = TRUE)
n_up  <- sum(results_df$padj < 0.05 & results_df$logFC > 0, na.rm = TRUE)
n_dn  <- sum(results_df$padj < 0.05 & results_df$logFC < 0, na.rm = TRUE)
cat("=== DEG summary (FDR < 0.05) ===\n")
cat("Significant:", n_sig, " | Up:", n_up, " | Down:", n_dn, "\n\n")

# ------------------ 11) Volcano plot ------------------
log2fc_cutoff <- 1
padj_cutoff <- 0.05

res_df <- results_df %>%
  filter(!is.na(padj) & !is.na(logFC)) %>%
  rename(log2FoldChange = logFC) %>%
  mutate(
    significance = case_when(
      log2FoldChange >  log2fc_cutoff & padj < padj_cutoff ~ "up",
      log2FoldChange < -log2fc_cutoff & padj < padj_cutoff ~ "down",
      TRUE ~ "ns"
    )
  )

# Label a few strongest hits (stringent threshold)
top_genes <- res_df %>%
  filter(padj < 1e-10 & abs(log2FoldChange) > 2) %>%
  slice_max(order_by = -log10(padj) * abs(log2FoldChange), n = 30)

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_text_repel(data = top_genes, aes(label = gene), size = 3, max.overlaps = 30) +
  scale_color_manual(values = c("up" = "red", "down" = "blue", "ns" = "grey")) +
  geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed") +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
  coord_cartesian(xlim = c(-7, 7)) +
  theme_minimal() +
  labs(title = "Volcano Plot: NASH vs Healthy",
       x = "log2 Fold Change (NB-GLM coefficient)",
       y = "-log10 adjusted p-value",
       color = "Regulation")

# ------------------ 12) Save results + Top DEGs table ------------------
# Base mean from normalized counts (for interpretability)
base_means <- rowMeans(norm_counts)

# ensure it's a tibble and has a 'gene' column
results_df <- tibble::as_tibble(results_df)
if (!"gene" %in% names(results_df) && !is.null(rownames(results_df))) {
  results_df <- dplyr::mutate(results_df, gene = rownames(results_df))
}

deg_table <- results_df %>%
  dplyr::filter(!is.na(logFC) & !is.na(pval) & !is.na(padj)) %>%
  dplyr::mutate(
    baseMean = base_means[gene],
    stat = logFC / se
  ) %>%
  dplyr::rename(
    log2FoldChange = logFC,
    lfcSE = se,
    pvalue = pval
  ) %>%
  dplyr::select(gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>%
  dplyr::arrange(padj)


# Write CSV
out_path <- "nb_glm_differential_expression_results.csv"
write.csv(deg_table, file = out_path, row.names = FALSE)
cat("Results written to:", out_path, "\n\n")

# Show top hits in console
cat("=== Top 10 DEGs (by FDR) ===\n")
print(head(deg_table, 10))
cat("\n")

# ------------------ 13) Integrity checks ------------------
stopifnot(all(names(sample_size_factors) == colnames(counts_filt)))
stopifnot(all.equal(rownames(norm_counts), rownames(counts_filt)))
