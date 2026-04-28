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
  library(scales)
})

# File paths (assumes running from project root)
all_reads_file <- "gene_raw_counts_all_reads.tsv"
hg_reads_file <- "gene_raw_counts_hg_reads.tsv"
output_dir <- "output_plots"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 1. Read the data
cat("Reading matrices...\n")
# read_tsv is generally faster and more tidyverse-idiomatic
df_all <- read_tsv(all_reads_file, show_col_types = FALSE)
df_hg <- read_tsv(hg_reads_file, show_col_types = FALSE)

# --- Filter Spike-ins ---
# Ignore records where gene_id starts with "ERCC-" or "SIRV"
cat("Filtering out ERCC and SIRV records...\n")
df_all <- df_all %>% filter(!grepl("^(ERCC-|SIRV)", gene_id))
df_hg <- df_hg %>% filter(!grepl("^(ERCC-|SIRV)", gene_id))
# ------------------------

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

# --- Verification Test on Input Data ---
# Ensure all common sample IDs contain at least one underscore for dataset extraction.
cat("Verifying sample ID formats...\n")
ids_without_underscore <- common_samples[!grepl("_", common_samples)]
if (length(ids_without_underscore) > 0) {
  stop("Dataset extraction regex test failed! The following common sample IDs do not contain an underscore:\n",
       paste(head(ids_without_underscore, 10), collapse=", "), 
       if(length(ids_without_underscore) > 10) " ..." else "")
}
# ---------------------------------------

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

# Process each sample and collect results
cat("Processing individual samples and generating plots...\n")
correlation_results <- map_dfr(common_samples, function(s_id) {
  sample_data <- df_combined %>%
    filter(sample_id == s_id) %>%
    filter(!is.na(counts_all), !is.na(counts_hg))
  
  if (nrow(sample_data) < 2) {
    cat("Skipping sample", s_id, ": insufficient data.\n")
    return(NULL)
  }
  
  # Calculate correlations on log-transformed counts
  log_all <- log10(sample_data$counts_all + 1)
  log_hg <- log10(sample_data$counts_hg + 1)
  
  r_pearson  <- cor(log_all, log_hg, method = "pearson")
  r_spearman <- cor(sample_data$counts_all, sample_data$counts_hg, method = "spearman")
  
  cat("Sample:", s_id, "| Pearson R (log10):", round(r_pearson, 4), "| Spearman Rho:", round(r_spearman, 4), "\n")
  
  # Create individual plot
  p <- ggplot(sample_data, aes(x = counts_all + 1, y = counts_hg + 1)) +
    geom_point(alpha = 0.2, size = 0.5) +
    scale_x_log10(labels = label_scientific()) +
    scale_y_log10(labels = label_scientific()) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    labs(
      title = paste("Gene Count Correlation -", s_id),
      subtitle = paste0("Pearson R (log10) = ", round(r_pearson, 4), 
                        "\nSpearman Rho = ", round(r_spearman, 4),
                        "\n(n = ", nrow(sample_data), " genes)"),
      x = "Raw Counts + 1 (All Reads, log10)",
      y = "Raw Counts + 1 (HG Reads, log10)"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  # Save as PNG
  file_name <- file.path(output_dir, paste0(s_id, "_correlation.png"))
  ggsave(file_name, plot = p, width = 7, height = 7, dpi = 150)
  
  # Return data for summary plot
  return(data.frame(
    sample_id = s_id,
    pearson_r = r_pearson,
    spearman_rho = r_spearman,
    stringsAsFactors = FALSE
  ))
})

# Generate Summary Boxplot
if (nrow(correlation_results) > 0) {
  cat("\nGenerating summary boxplot...\n")
  
  # Extract dataset name (everything until the last underscore)
  # All IDs guaranteed to have an underscore by the test at start
  correlation_results <- correlation_results %>%
    mutate(dataset = sub("_[^_]*$", "", sample_id))
  
  p_summary <- ggplot(correlation_results, aes(x = dataset, y = pearson_r, fill = dataset)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
    labs(
      title = "Pearson Correlation Summary by Dataset",
      subtitle = "Correlations calculated on log10(counts + 1)",
      x = "Dataset",
      y = "Pearson R (log10 counts)"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  summary_file <- file.path(output_dir, "dataset_pearson_summary.png")
  ggsave(summary_file, plot = p_summary, width = 10, height = 7, dpi = 150)
  cat("Summary plot saved to:", summary_file, "\n")
}

cat("\nAll plots have been saved to:", output_dir, "\n")
cat("Done.\n")
