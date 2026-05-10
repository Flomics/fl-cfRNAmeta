#!/usr/bin/env Rscript

# src/correlation_analysis.R
# Processes gene count matrices sample by sample
# Produces individual scatterplots and summary plots of correlations.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(purrr)
  library(scales)
  library(jsonlite)
})

# --- Robust Data Loading Helper ---
robust_read <- function(file, name, n_max = Inf) {
  df <- read_tsv(file, show_col_types = FALSE, guess_max = 100000, quote = "", n_max = n_max)
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

mappings_json <- "fl-cfRNAmeta/src/dataset_mappings.json"

if (!dir.exists(output_dir)) {
  cat("Creating output directory:", output_dir, "\n")
  dir.create(output_dir, recursive = TRUE)
}

# Metadata
df_mapping <- robust_read(mapping_file, "mapping file")
df_sampleinfo <- robust_read(sampleinfo_file, "sampleinfo file") %>%
  select(sample_name, avg_mapped_read_length, mapped_percentage)

if (!only_summary_plot) {
  # -----------------------------------------
  # --- Full Analysis Mode ---
  # -----------------------------------------
  cat("Reading matrices...\n")
  df_all <- robust_read(all_reads_file, "All Reads matrix")
  df_hg <- robust_read(hg_reads_file, "HG Reads matrix")

  df_all <- df_all %>% filter(!grepl("^(ERCC-|SIRV)", gene_id))
  df_hg <- df_hg %>% filter(!grepl("^(ERCC-|SIRV)", gene_id))

  all_header <- robust_read(all_reads_file, "All Reads header", n_max = 0)
  hg_header  <- robust_read(hg_reads_file, "HG Reads header", n_max = 0)
  common_samples <- intersect(colnames(all_header)[-(1:2)], colnames(hg_header)[-(1:2)])

  if (!is.null(samples_to_keep_file)) {
    target_samples <- read_lines(samples_to_keep_file) %>% trimws() %>% .[. != ""]
    common_samples <- intersect(common_samples, target_samples)
  }

  cat("Processing", length(common_samples), "common samples...\n")

  df_all_long <- df_all %>% select(gene_id, gene_name, all_of(common_samples)) %>%
    pivot_longer(cols = -c(gene_id, gene_name), names_to = "sample_id", values_to = "counts_all")
  df_hg_long <- df_hg %>% select(gene_id, gene_name, all_of(common_samples)) %>%
    pivot_longer(cols = -c(gene_id, gene_name), names_to = "sample_id", values_to = "counts_hg")
  df_combined <- inner_join(df_all_long, df_hg_long, by = c("gene_id", "gene_name", "sample_id"))

  correlation_results <- map_dfr(common_samples, function(s_id) {
    sample_data <- df_combined %>% filter(sample_id == s_id) %>% filter(!is.na(counts_all), !is.na(counts_hg))
    if (nrow(sample_data) < 2) return(NULL)
    r_pearson  <- cor(log10(sample_data$counts_all + 1), log10(sample_data$counts_hg + 1), method = "pearson")
    r_spearman <- cor(sample_data$counts_all, sample_data$counts_hg, method = "spearman")
    if (r_pearson < 0.9) {
      outliers <- sample_data %>% mutate(log_diff = abs(log10(counts_all + 1) - log10(counts_hg + 1))) %>%
        arrange(desc(log_diff)) %>% head(1000)
      write_tsv(outliers, file.path(output_dir, paste0(s_id, "_top1000outliers.tsv")))
    }
    # Individual plot
    p <- ggplot(sample_data, aes(x = counts_all + 1, y = counts_hg + 1)) +
      geom_point(alpha = 0.2, size = 0.5) + scale_x_log10(labels = label_scientific()) +
      scale_y_log10(labels = label_scientific()) + geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
      labs(title = paste("Correlation -", s_id), subtitle = paste0("R = ", round(r_pearson, 4)),
           x = "Counts+1 (All)", y = "Counts+1 (HG)") + theme_minimal()
    ggsave(file.path(output_dir, paste0(s_id, "_correlation.png")), plot = p, width = 7, height = 7, dpi = 150)
    return(data.frame(sample_id = s_id, pearson_r = r_pearson, spearman_rho = r_spearman, stringsAsFactors = FALSE))
  })

  correlation_results <- correlation_results %>% left_join(df_sampleinfo, by = c("sample_id" = "sample_name"))
  write_tsv(correlation_results, correlations_file)

} else {
  # -----------------------------------------
  # --- Summary Only Mode ---
  # -----------------------------------------
  correlation_results <- robust_read(correlations_file, "pre-calculated correlations")
  all_header <- robust_read(all_reads_file, "All Reads header", n_max = 0)
  hg_header  <- robust_read(hg_reads_file, "HG Reads header", n_max = 0)
  target_common <- intersect(colnames(all_header)[-(1:2)], colnames(hg_header)[-(1:2)])
  if (!is.null(samples_to_keep_file)) {
    targets <- read_lines(samples_to_keep_file) %>% trimws() %>% .[. != ""]
    target_common <- intersect(target_common, targets)
  }
  correlation_results <- correlation_results %>% filter(sample_id %in% target_common)
  if (!all(c("avg_mapped_read_length", "mapped_percentage") %in% colnames(correlation_results))) {
    correlation_results <- correlation_results %>% select(-any_of(c("avg_mapped_read_length", "mapped_percentage"))) %>%
      left_join(df_sampleinfo, by = c("sample_id" = "sample_name"))
  }
}

# -----------------------------------------
# --- Generate Plots ---
# -----------------------------------------
if (nrow(correlation_results) > 0) {
  plot_data <- correlation_results %>% left_join(df_mapping, by = c("sample_id" = "sample_name")) %>%
    rename(dataset = dataset_batch) %>% filter(!is.na(dataset), dataset != "")

  final_palette <- NULL
  if (file.exists(mappings_json)) {
    m_json <- fromJSON(mappings_json)
    v_order <- m_json$datasetVisualOrder[m_json$datasetVisualOrder %in% unique(plot_data$dataset)]
    final_order <- c(v_order, setdiff(unique(plot_data$dataset), v_order))
    v_labels <- unlist(m_json$datasetsLabels)
    if (!is.null(m_json$datasetsPalette)) {
      final_palette <- unlist(m_json$datasetsPalette)[names(unlist(m_json$datasetsPalette)) %in% names(v_labels)]
      names(final_palette) <- v_labels[names(final_palette)]
    }
    plot_data <- plot_data %>% mutate(dataset = factor(dataset, levels = final_order)) %>%
      mutate(dataset_label = ifelse(dataset %in% names(v_labels), v_labels[as.character(dataset)], as.character(dataset))) %>%
      mutate(dataset_label = factor(dataset_label, levels = v_labels[as.character(final_order)])) %>%
      mutate(dataset = dataset_label)
  }

  # 1. Summary Pearson Boxplot
  p_summary <- ggplot(plot_data, aes(x = dataset, y = pearson_r)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, fill = NA, color = "lightgrey") +
    geom_jitter(aes(color = avg_mapped_read_length), width = 0.2, alpha = 0.5, size = 1.5) +
    scale_y_continuous(limits = c(0, 1)) + scale_color_viridis_c(option = "viridis") +
    labs(title = "Pearson Correlation Summary by Dataset", x = "Dataset", y = "Pearson R (log10)", color = "Read Length") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
  ggsave(file.path(output_dir, "dataset_pearson_summary.png"), plot = p_summary, width = 14, height = 7, dpi = 150)

  # 2. Pearson vs Mapped Percentage (Faceted)
  facet_correlations <- plot_data %>% group_by(dataset) %>%
    summarize(r_val = cor(mapped_percentage, pearson_r, use = "complete.obs"), n_samples = n(), .groups = "drop") %>%
    mutate(label = paste0("r = ", round(r_val, 3), "\nn = ", n_samples))
  p_faceted <- ggplot(plot_data, aes(x = mapped_percentage, y = pearson_r)) +
    geom_point(aes(color = avg_mapped_read_length), alpha = 0.7, size = 2) +
    geom_text(data = facet_correlations, aes(x = Inf, y = 0, label = label), hjust = 1.1, vjust = -0.5, size = 3, inherit.aes = FALSE) +
    scale_y_continuous(limits = c(0, 1)) + scale_color_viridis_c(option = "viridis") +
    facet_wrap(~dataset, ncol = 6) + labs(title = "Pearson R vs Mapped % by Dataset", x = "Mapped %", y = "Pearson R (log10)", color = "Read Length") +
    theme_minimal() + theme(legend.position = "bottom")
  ggsave(file.path(output_dir, "all_datasets_mapped_pct_vs_pearson.png"), plot = p_faceted, width = 18, height = 3 * ceiling(length(unique(plot_data$dataset))/6) + 2, dpi = 150)

  # 3. Pearson vs Avg Mapped Read Length (Global) + STATS
  cat("\nGenerating Pearson vs Avg Mapped Read Length plot and statistical analysis...\n")
  global_r <- cor(plot_data$avg_mapped_read_length, plot_data$pearson_r, use = "complete.obs")
  p_rl <- ggplot(plot_data, aes(x = avg_mapped_read_length, y = pearson_r)) +
    geom_point(aes(color = dataset), alpha = 0.6, size = 2) +
    geom_smooth(method = "loess", color = "black", se = FALSE, linetype = "solid", linewidth = 0.8) +
    scale_y_continuous(limits = c(0, 1)) + labs(title = "", x = "Effective fragment length\n(average mapped length, bp)", y = "Pearson R", color = "Dataset") +
    theme_minimal() + theme(legend.position = "right")
  if (!is.null(final_palette)) p_rl <- p_rl + scale_color_manual(values = final_palette)
  ggsave(file.path(output_dir, "pearson_vs_read_length.png"), plot = p_rl, width = 12, height = 7, dpi = 150)

  # --- Threshold Statistical Analysis (X = 100) ---
  threshold <- 100
  stat_data <- plot_data %>% filter(!is.na(avg_mapped_read_length), !is.na(pearson_r)) %>%
    mutate(group = ifelse(avg_mapped_read_length < threshold, paste0("< ", threshold), paste0(">= ", threshold)))
  
  group_summary <- stat_data %>% group_by(group) %>%
    summarize(n = n(), mean_r = mean(pearson_r), median_r = median(pearson_r), sd_r = sd(pearson_r), .groups = "drop")
  
  # Wilcoxon rank sum test (non-parametric)
  wilcox_res <- wilcox.test(pearson_r ~ group, data = stat_data)
  
  stats_file <- file.path(output_dir, "read_length_threshold_analysis.txt")
  sink(stats_file)
  cat("========================================================================\n")
  cat("Statistical Analysis: Pearson R by Read Length Threshold (", threshold, "bp)\n", sep="")
  cat("========================================================================\n\n")
  cat("Group Summary Statistics:\n")
  print(as.data.frame(group_summary))
  cat("\n------------------------------------------------------------------------\n")
  cat("Hypothesis: Pearson R values differ between the two groups.\n\n")
  cat("1. Wilcoxon Rank Sum Test with Continuity Correction (Non-Parametric):\n")
  print(wilcox_res)
  sink()
  cat("Statistical analysis results saved to:", stats_file, "\n")
}

cat("\nDone.\n")
