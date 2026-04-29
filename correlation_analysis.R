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

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4 || length(args) > 5) {
  cat("Usage: correlation_analysis.R <all_reads_file> <hg_reads_file> <mapping_file> <output_dir> [samples_to_keep_file]\n")
  quit(save = "no", status = 1)
}

all_reads_file <- args[1]
hg_reads_file  <- args[2]
mapping_file   <- args[3]
output_dir     <- args[4]
samples_to_keep_file <- if (length(args) == 5) args[5] else NULL
correlations_file <- file.path(output_dir, "correlations.tsv")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 1. Read the data
cat("Reading matrices...\n")
# read_tsv is generally faster and more tidyverse-idiomatic
df_all <- read_tsv(all_reads_file, show_col_types = FALSE)
df_hg <- read_tsv(hg_reads_file, show_col_types = FALSE)
df_mapping <- read_tsv(mapping_file, show_col_types = FALSE)

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

# --- Filter Samples by List (Optional) ---
if (!is.null(samples_to_keep_file)) {
  cat("Filtering samples using list from:", samples_to_keep_file, "\n")
  if (!file.exists(samples_to_keep_file)) {
    cat("ERROR: Sample list file not found:", samples_to_keep_file, "\n")
    quit(save = "no", status = 1)
  }
  # Read the list of samples, stripping whitespace and skipping empty lines
  target_samples <- read_lines(samples_to_keep_file) %>% trimws() %>% .[. != ""]
  
  original_n <- length(common_samples)
  common_samples <- intersect(common_samples, target_samples)
  
  if (length(common_samples) == 0) {
    cat("ERROR: No common samples remain after filtering with the provided list.\n")
    quit(save = "no", status = 1)
  }
  cat("Kept", length(common_samples), "out of", original_n, "samples based on the provided list.\n")
}
# -----------------------------------------

# --- Check Mapping Presence ---
cat("Checking if all samples have mapping information...\n")
missing_metadata <- setdiff(common_samples, df_mapping$sample_name)
if (length(missing_metadata) > 0) {
  cat("ERROR: The following samples are missing from the mapping file:", mapping_file, "\n")
  cat(paste(missing_metadata, collapse = ", "), "\n")
  quit(save = "no", status = 1)
}
# ------------------------------

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

# Save correlations to TSV
cat("Saving correlation results to:", correlations_file, "\n")
write_tsv(correlation_results, correlations_file)

# Generate Summary Boxplot
if (nrow(correlation_results) > 0) {
  cat("\nGenerating summary boxplot...\n")
  
  # Map sample names to dataset names using df_mapping
  correlation_results <- correlation_results %>%
    left_join(df_mapping, by = c("sample_id" = "sample_name")) %>%
    rename(dataset = dataset_batch)
  
  # Identify and warn about samples with empty dataset information
  samples_with_na_dataset <- correlation_results %>%
    filter(is.na(dataset) | dataset == "")
  
  if (nrow(samples_with_na_dataset) > 0) {
    cat("WARNING: The following samples have no dataset mapping and will be excluded from the summary plot:\n")
    cat(paste(samples_with_na_dataset$sample_id, collapse = ", "), "\n")
    correlation_results <- correlation_results %>%
      filter(!is.na(dataset), dataset != "")
  }

  if (nrow(correlation_results) > 0) {
    p_summary <- ggplot(correlation_results, aes(x = dataset, y = pearson_r, fill = dataset)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
      scale_y_continuous(limits = c(0, 1)) +
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
  } else {
    cat("WARNING: No samples with valid dataset mapping remaining. Skipping summary plot.\n")
  }
}

cat("\nAll plots have been saved to:", output_dir, "\n")
cat("Done.\n")
