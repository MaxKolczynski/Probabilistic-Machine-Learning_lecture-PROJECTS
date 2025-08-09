
#install.packages("BiocManager")
#BiocManager::install("DESeq2")

# Load libraries
suppressMessages(library("DESeq2"))



# Load data
count_data <- read.csv("pseudobulk_counts.csv", row.names = 1)
meta_data <- read.csv("sample_metadata.csv", row.names = 1)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = meta_data,
  design = ~ condition
)

# Run DE analysis
dds <- DESeq(dds)
res <- results(dds)

# Order by adjusted p-value
resOrdered <- res[order(res$padj), ]

# Save result
write.csv(as.data.frame(resOrdered), file = "deseq2_results.csv")