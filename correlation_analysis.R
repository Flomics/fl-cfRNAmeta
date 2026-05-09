#!/usr/bin/env Rscript

# src/correlation_analysis.R
# Processes gene count matrices sample by sample
# Produces individual scatterplots and a summary boxplot of correlations.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(purrr)
  library(scales)
})

# --- Robust Data Loading Helper ---
robust_read <- function(file, name, n_max = Inf) {
  # Increase guess_max for large files and disable quoting
  df <- read_tsv(file, show_col_types = FALSE, guess_max = 100000, quote = "", n_max = n_max)
  
  # Report parsing problems if any
  p <- problems(df)
  if (nrow(p) > 0) {
    cat(paste0("\nWARNING: Parsing issues detected in ", name, " (", file, "):\n"), file = stderr())
    print(head(p, 5), file = stderr())
    cat("... use problems(df) in R for the full list.\n\n", file = stderr())
  }
  return(df)
}

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Parse flags
only_summary_plot <- "--only_summary_plot" %in% args
args <- args[args != "--only_summary_plot"]

if (length(args) < 5 || length(args) > 6) {
  cat("Usage: correlation_analysis.R [--only_summary_plot] <all_reads_file> <hg_reads_file> <mapping_file> <sampleinfo_file> <output_dir> [samples_to_keep_file]\n")
  quit(save = "no", status = 1)
}

all_reads_file    <- args[1]
hg_reads_file     <- args[2]
mapping_file      <- args[3]
sampleinfo_file   <- args[4]
output_dir        <- args[5]
samples_to_keep_file <- if (length(args) == 6) args[6] else NULL
correlations_file <- file.path(output_dir, "correlations.tsv")

if (!dir.exists(output_dir)) {
  cat("Creating output directory:", output_dir, "\n")
  dir.create(output_dir, recursive = TRUE)
}

# --- 1. Identify Expected Samples (Efficiency: Read only headers) ---
cat("Identifying expected samples from matrix headers...\n")
all_header <- robust_read(all_reads_file, "All Reads header", n_max = 0)
hg_header  <- robust_read(hg_reads_file, "HG Reads header", n_max = 0)

all_samples <- colnames(all_header)[-(1:2)]
hg_samples  <- colnames(hg_header)[-(1:2)]
common_samples <- intersect(all_samples, hg_samples)

if (length(common_samples) == 0) {
  cat("\nERROR: No common sample IDs found between the two matrices.\n")
  quit(save = "no", status = 1)
}

# --- 2. Filter Samples by List (Optional) ---
if (!is.null(samples_to_keep_file)) {
  cat("Filtering samples using list from:", samples_to_keep_file, "\n")
  if (!file.exists(samples_to_keep_file)) {
    cat("ERROR: Sample list file not found:", samples_to_keep_file, "\n")
    quit(save = "no", status = 1)
  }
  target_samples <- read_lines(samples_to_keep_file) %>% trimws() %>% .[. != ""]
  original_n <- length(common_samples)
  common_samples <- intersect(common_samples, target_samples)
  if (length(common_samples) == 0) {
    cat("ERROR: No common samples remain after filtering.\n")
    quit(save = "no", status = 1)
  }
  cat("Kept", length(common_samples), "out of", original_n, "samples based on provided list.\n")
}

# --- 3. Read Metadata ---
df_mapping <- robust_read(mapping_file, "mapping file")
cat("Reading sampleinfo from:", sampleinfo_file, "\n")
df_sampleinfo <- robust_read(sampleinfo_file, "sampleinfo file") %>%
  select(sample_name, avg_mapped_read_length, mapped_percentage)

# Check Mapping Presence
missing_metadata <- setdiff(common_samples, df_mapping$sample_name)
if (length(missing_metadata) > 0) {
  cat("ERROR: The following samples are missing from mapping file:\n")
  cat(paste(missing_metadata, collapse = ", "), "\n")
  quit(save = "no", status = 1)
}


if (!only_summary_plot) {
  # -----------------------------------------
  # --- Full Analysis Mode ---
  # -----------------------------------------
  cat("Reading matrices...\n")
  df_all <- robust_read(all_reads_file, "All Reads matrix")
  df_hg <- robust_read(hg_reads_file, "HG Reads matrix")

  # --- Filter Spike-ins ---
  cat("Filtering out ERCC and SIRV records...\n")
  df_all <- df_all %>% filter(!grepl("^(ERCC-|SIRV)", gene_id))
  df_hg <- df_hg %>% filter(!grepl("^(ERCC-|SIRV)", gene_id))

  cat("Processing", length(common_samples), "common samples...\n")

  cat("Reshaping data...\n")
  df_all_long <- df_all %>%
    select(gene_id, gene_name, all_of(common_samples)) %>%
    pivot_longer(cols = -c(gene_id, gene_name), names_to = "sample_id", values_to = "counts_all")

  df_hg_long <- df_hg %>%
    select(gene_id, gene_name, all_of(common_samples)) %>%
    pivot_longer(cols = -c(gene_id, gene_name), names_to = "sample_id", values_to = "counts_hg")

  df_combined <- inner_join(df_all_long, df_hg_long, by = c("gene_id", "gene_name", "sample_id"))

  cat("Calculating correlations and generating plots...\n")
  correlation_results <- map_dfr(common_samples, function(s_id) {
    sample_data <- df_combined %>%
      filter(sample_id == s_id) %>%
      filter(!is.na(counts_all), !is.na(counts_hg))
    
    if (nrow(sample_data) < 2) return(NULL)
    
    log_all <- log10(sample_data$counts_all + 1)
    log_hg <- log10(sample_data$counts_hg + 1)
    
    r_pearson  <- cor(log_all, log_hg, method = "pearson")
    r_spearman <- cor(sample_data$counts_all, sample_data$counts_hg, method = "spearman")
    
    if (r_pearson < 0.9) {
      cat("  Sample", s_id, "has low correlation (R =", round(r_pearson, 4), "). Saving top 1000 outliers...\n")
      outliers <- sample_data %>%
        mutate(log_diff = abs(log10(counts_all + 1) - log10(counts_hg + 1))) %>%
        arrange(desc(log_diff)) %>%
        head(1000)
      write_tsv(outliers, file.path(output_dir, paste0(s_id, "_top1000outliers.tsv")))
    }

    # Individual plot
    p <- ggplot(sample_data, aes(x = counts_all + 1, y = counts_hg + 1)) +
      geom_point(alpha = 0.2, size = 0.5) +
      scale_x_log10(labels = label_scientific()) +
      scale_y_log10(labels = label_scientific()) +
      geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
      labs(
        title = paste("Gene Count Correlation -", s_id),
        subtitle = paste0("Pearson R (log10) = ", round(r_pearson, 4), 
                          "\nSpearman Rho = ", round(r_spearman, 4)),
        x = "Raw Counts + 1 (All Reads, log10)",
        y = "Raw Counts + 1 (HG Reads, log10)"
      ) + theme_minimal()
    
    ggsave(file.path(output_dir, paste0(s_id, "_correlation.png")), plot = p, width = 7, height = 7, dpi = 150)
    
    return(data.frame(sample_id = s_id, pearson_r = r_pearson, spearman_rho = r_spearman, stringsAsFactors = FALSE))
  })

  cat("Saving correlation results to:", correlations_file, "\n")
  write_tsv(correlation_results, correlations_file)

} else {
  # -----------------------------------------
  # --- Summary Only Mode ---
  # -----------------------------------------
  cat("Mode: --only_summary_plot. Reading pre-calculated results from:", correlations_file, "\n")
  if (!file.exists(correlations_file)) {
    cat("ERROR: Correlations file not found. Run without --only_summary_plot first.\n")
    quit(save = "no", status = 1)
  }
  correlation_results <- robust_read(correlations_file, "pre-calculated correlations")

  # Ensure all common_samples are present in loaded correlations
  missing_correlations <- setdiff(common_samples, correlation_results$sample_id)
  if (length(missing_correlations) > 0) {
    cat("ERROR: Pre-calculated correlations are missing for the following samples:\n")
    cat(paste(missing_correlations, collapse = ", "), "\n")
    cat("Please run the full analysis without --only_summary_plot to generate them.\n")
    quit(save = "no", status = 1)
  }
  
  # Filter the loaded results to only include the common_samples (in case TSV has more)
  correlation_results <- correlation_results %>% filter(sample_id %in% common_samples)
  cat("Verified and kept", nrow(correlation_results), "samples for summary plotting.\n")
}

# -----------------------------------------
# --- Generate Plots ---
# -----------------------------------------
if (nrow(correlation_results) > 0) {
  
  # Join all metadata for plotting
  plot_data <- correlation_results %>%
    left_join(df_mapping, by = c("sample_id" = "sample_name")) %>%
    rename(dataset = dataset_batch) %>%
    filter(!is.na(dataset), dataset != "") %>%
    left_join(df_sampleinfo, by = c("sample_id" = "sample_name"))

  if (nrow(plot_data) > 0) {
    
    # 1. Summary Pearson Boxplot
    cat("\nGenerating summary boxplot...\n")
    p_summary <- ggplot(plot_data, aes(x = dataset, y = pearson_r)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA, fill = NA, color = "lightgrey") +
      geom_jitter(aes(color = avg_mapped_read_length), width = 0.2, alpha = 0.5, size = 1.5) +
      scale_y_continuous(limits = c(0, 1)) +
      scale_color_viridis_c(option = "viridis") +
      labs(
        title = "Pearson Correlation Summary by Dataset",
        subtitle = "Correlations calculated on log10(counts + 1)",
        x = "Dataset",
        y = "Pearson R (log10 counts)",
        color = "Avg Mapped Read Length"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "right"
      )
    
    summary_file <- file.path(output_dir, "dataset_pearson_summary.png")
    ggsave(summary_file, plot = p_summary, width = 12, height = 7, dpi = 150)
    cat("Summary plot saved to:", summary_file, "\n")

    # 2. Pearson vs Mapped Percentage (One per Dataset)
    cat("\nGenerating Pearson vs Mapped Percentage scatterplots per dataset...\n")
    datasets <- unique(plot_data$dataset)
    for (ds in datasets) {
      ds_data <- plot_data %>% filter(dataset == ds)
      
      p_ds <- ggplot(ds_data, aes(x = mapped_percentage, y = pearson_r)) +
        geom_point(aes(color = avg_mapped_read_length), alpha = 0.7, size = 3) +
        scale_y_continuous(limits = c(0, 1)) +
        scale_color_viridis_c(option = "viridis") +
        labs(
          title = paste("Pearson R vs Mapped % -", ds),
          subtitle = "Correlations calculated on log10(counts + 1)",
          x = "Mapped Percentage (%)",
          y = "Pearson R (log10 counts)",
          color = "Avg Mapped Read Length"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "right"
        )
      
      ds_plot_file <- file.path(output_dir, paste0(ds, "_mapped_pct_vs_pearson.png"))
      ggsave(ds_plot_file, plot = p_ds, width = 8, height = 7, dpi = 150)
      cat("  Saved dataset plot:", ds_plot_file, "\n")
    }

  } else {
    cat("WARNING: No samples with valid dataset mapping remaining. Skipping plots.\n")
  }
}

cat("\nDone.\n")
