#!/usr/bin/env Rscript

# src/correlation_analysis.R
# Processes gene count matrices sample by sample using tidyverse.
# Produces one scatterplot per sample with Pearson correlation.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(purrr)
})

# File paths (assumes running from project root)
all_reads_file <- "gene_raw_counts_all_reads.tsv"
hg_reads_file <- "gene_raw_counts_hg_reads.tsv"
output_dir <- "output_plots"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cat("Reading matrices...\n")
# read_tsv is generally faster and more tidyverse-idiomatic
df_all <- read_tsv(all_reads_file, show_col_types = FALSE)
df_hg <- read_tsv(hg_reads_file, show_col_types = FALSE)

# Identify common samples (columns 3 onwards)
all_samples <- colnames(df_all)[-(1:2)]
hg_samples <- colnames(df_hg)[-(1:2)]
common_samples <- intersect(all_samples, hg_samples)

if (length(common_samples) == 0) {
  cat("\nWARNING: No common sample IDs found between the two files.\n")
  cat("Samples in All Reads (first 3): ", paste(head(all_samples, 3), collapse=", "), "\n")
  cat("Samples in HG Reads (first 3): ", paste(head(hg_samples, 3), collapse=", "), "\n")
  quit(save = "no", status = 1)
}

cat("Processing", length(common_samples), "common samples...\n")

# Reshape to long format for easier joining and per-sample processing
cat("Reshaping data...\n")
df_all_long <- df_all %>%
  select(gene_id, all_of(common_samples)) %>%
  pivot_longer(cols = -gene_id, names_to = "sample_id", values_to = "counts_all")

df_hg_long <- df_hg %>%
  select(gene_id, all_of(common_samples)) %>%
  pivot_longer(cols = -gene_id, names_to = "sample_id", values_to = "counts_hg")

# Join the two datasets
df_combined <- inner_join(df_all_long, df_hg_long, by = c("gene_id", "sample_id"))

# Process each sample
walk(common_samples, function(s_id) {
  sample_data <- df_combined %>%
    filter(sample_id == s_id) %>%
    filter(!is.na(counts_all), !is.na(counts_hg))
  
  if (nrow(sample_data) < 2) {
    cat("Skipping sample", s_id, ": insufficient data.\n")
    return()
  }
  
  # Calculate Pearson correlation
  r_val <- cor(sample_data$counts_all, sample_data$counts_hg, method = "pearson")
  
  cat("Sample:", s_id, "| Pearson R:", round(r_val, 4), "\n")
  
  # Create plot
  p <- ggplot(sample_data, aes(x = counts_all + 1, y = counts_hg + 1)) +
    geom_point(alpha = 0.2, size = 0.5) +
    scale_x_log10() +
    scale_y_log10() +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    labs(
      title = paste("Gene Count Correlation -", s_id),
      subtitle = paste0("Pearson R = ", round(r_val, 4), " (n = ", nrow(sample_data), " genes)"),
      x = "Raw Counts + 1 (All Reads, log10)",
      y = "Raw Counts + 1 (HG Reads, log10)"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  # Save as PNG
  file_name <- file.path(output_dir, paste0(s_id, "_correlation.png"))
  ggsave(file_name, plot = p, width = 7, height = 7, dpi = 150)
})

cat("\nAll plots have been saved to:", output_dir, "\n")
cat("Done.\n")
