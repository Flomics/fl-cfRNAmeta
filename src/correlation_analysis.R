#!/usr/bin/env Rscript

# correlation_analysis.R
# Script to process gene count matrices and produce a correlation scatterplot.

# File paths
all_reads_file <- "gene_raw_counts_all_reads.tsv"
hg_reads_file <- "gene_raw_counts_hg_reads.tsv"
output_plot <- "gene_count_correlation.pdf"

# 1. Read the data
cat("Reading matrices...\n")
# Using check.names=FALSE to preserve sample IDs (like SRR... or block_...)
df_all <- read.delim(all_reads_file, check.names = FALSE, stringsAsFactors = FALSE)
df_hg <- read.delim(hg_reads_file, check.names = FALSE, stringsAsFactors = FALSE)

# 2. Identify common samples and genes
# Samples are columns 3 to N
all_samples <- colnames(df_all)[-(1:2)]
hg_samples <- colnames(df_hg)[-(1:2)]
common_samples <- intersect(all_samples, hg_samples)

if (length(common_samples) == 0) {
  cat("\nWARNING: No common sample IDs found between the two files.\n")
  cat("Common columns found: ", paste(intersect(colnames(df_all), colnames(df_hg)), collapse=", "), "\n")
  cat("Samples in All Reads (first 3): ", paste(head(all_samples, 3), collapse=", "), "\n")
  cat("Samples in HG Reads (first 3): ", paste(head(hg_samples, 3), collapse=", "), "\n")
  cat("\nCheck if sample IDs match or if a mapping is required.\n")
  quit(save = "no", status = 1)
}

cat("Found", length(common_samples), "common samples.\n")

# Genes (assuming gene_id is unique)
common_genes <- intersect(df_all$gene_id, df_hg$gene_id)
if (length(common_genes) == 0) {
  cat("ERROR: No common gene IDs found.\n")
  quit(save = "no", status = 1)
}
cat("Found", length(common_genes), "common genes.\n")

# 3. Align and subset matrices
# Use gene_id as row names for easy subsetting
rownames(df_all) <- df_all$gene_id
rownames(df_hg) <- df_hg$gene_id

# Extract matching data
mat_all <- as.matrix(df_all[common_genes, common_samples])
mat_hg <- as.matrix(df_hg[common_genes, common_samples])

# 4. Calculate correlation
# Flatten matrices to vectors
vals_x <- as.vector(mat_all)
vals_y <- as.vector(mat_hg)

# Remove any potential NAs
valid_idx <- !is.na(vals_x) & !is.na(vals_y)
vals_x <- vals_x[valid_idx]
vals_y <- vals_y[valid_idx]

pearson_corr <- cor(vals_x, vals_y, method = "pearson")
cat("Pearson Correlation:", round(pearson_corr, 4), "\n")

# 5. Generate plot
cat("Generating plot: ", output_plot, "...\n")
pdf(output_plot, width = 8, height = 8)

# Set up the plot with some transparency if many points
plot(vals_x, vals_y, 
     pch = 16, 
     col = rgb(0, 0, 0, 0.1), # Black with alpha
     main = "Gene Raw Counts Correlation",
     xlab = "Raw Counts (All Reads)",
     ylab = "Raw Counts (HG Reads)",
     cex = 0.5)

# Add identity line
abline(0, 1, col = "red", lty = 2)

# Add correlation text and info
legend("topleft", 
       legend = c(paste0("Pearson R = ", round(pearson_corr, 4)),
                  paste0("n_samples = ", length(common_samples)),
                  paste0("n_genes = ", length(common_genes))),
       bty = "n", cex = 1.2, text.col = "blue")

dev.off()
cat("Done.\n")
