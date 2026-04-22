library(dplyr)
library(tidyr)

######################################################
# Variance Partition Analysis
# Needs to run after the boxplots_fig2.R script
######################################################

suppressMessages(library("variancePartition"))
suppressMessages(library("edgeR"))


platelet_info <- read.delim("tables/sampleinfo_external_and_internal_datasets.tsv", header = TRUE, sep = "\t", fileEncoding = "UTF-8")

# Add platelet info to filtered_df so it's available for vp_sampleinfo
filtered_df <- filtered_df %>%
  left_join(
    platelet_info %>% select(sample_name, platelet),
    by = "sample_name"
  )

# ─── Datasets to include in the analysis ─────────────────────────────────────
DATASETS_TO_INCLUDE <- c(
  # Fill with values from dataset_batch.y to include
)

filtered_df <- filtered_df %>%
  filter(dataset_batch.y %in% DATASETS_TO_INCLUDE)

cat("Datasets included:", paste(DATASETS_TO_INCLUDE, collapse = ", "), "\n")
cat("Samples after dataset filter:", nrow(filtered_df), "\n")

# ─── Variables to test ───────────────────────────────────────────────────────
VP_NUMERIC <- c(
  "genes_contributing_to_80._of_reads",   # NG80
  "percentage_of_spliced_reads",           # FSR
  "exonic_reads_minus_spike_ins"           # FER 
)

VP_CATEGORICAL <- c(
  "dataset_batch.y",      # dataset
  "status"                # phenotype 
)

VP_NUMERIC <- c(
  "genes_contributing_to_80._of_reads",   # NG80
  "percentage_of_spliced_reads",           # FSR
  "exonic_reads_minus_spike_ins",          # FER
  "protein_coding_pct",                     # biotype composition
  "platelet"
)

VP_CATEGORICAL <- c(
  "status"            # phenotype
)

# ─── Prepare sampleinfo ──────────────────────────────────────────────────────
# Start from table_filtered which already has all the QC metrics
vp_sampleinfo <- table_filtered %>%
  select(sample_name, all_of(VP_NUMERIC), all_of(VP_CATEGORICAL)) %>%
  filter(!is.na(genes_contributing_to_80._of_reads) &
           !is.na(percentage_of_spliced_reads) &
           !is.na(status)) %>%
  as.data.frame()

vp_sampleinfo <- filtered_df %>%
  select(sample_name, dataset_batch.y, all_of(VP_NUMERIC), all_of(VP_CATEGORICAL)) %>%
  # Only require the complete variables — leave NAs for borderline ones
  filter(
    !is.na(genes_contributing_to_80._of_reads),
    !is.na(percentage_of_spliced_reads),
    !is.na(status)
  ) %>%
  as.data.frame()

# Report how many samples remain
cat("Samples after NA filtering:", nrow(vp_sampleinfo), "\n")
cat("Samples dropped:", nrow(filtered_df) - nrow(vp_sampleinfo), "\n")

row.names(vp_sampleinfo) <- vp_sampleinfo$sample_name

# Make categorical variables factors
for (cat in VP_CATEGORICAL) {
  vp_sampleinfo[[cat]] <- as.factor(vp_sampleinfo[[cat]])
}

# Remove categories with only 1 level (would break the model)
to_remove_vars <- c()
for (cat in VP_CATEGORICAL) {
  n <- length(unique(vp_sampleinfo[[cat]][!is.na(vp_sampleinfo[[cat]])]))
  if (n == 1) {
    message("Removing ", cat, " — only 1 unique level")
    to_remove_vars <- c(to_remove_vars, cat)
  }
}
VP_CATEGORICAL <- VP_CATEGORICAL[!VP_CATEGORICAL %in% to_remove_vars]
vp_sampleinfo  <- vp_sampleinfo[, !colnames(vp_sampleinfo) %in% to_remove_vars]

# Scale numeric variables (important for variance partition)
for (num in VP_NUMERIC) {
  vp_sampleinfo[[num]] <- scale(vp_sampleinfo[[num]])
}

# ─── Check collinearity between metadata variables ────────────────────────────
form_check <- as.formula(
  paste0(
    "~ ",
    paste(VP_NUMERIC, collapse = " + "),
    " + (1 | ", paste(VP_CATEGORICAL, collapse = ") + (1 | "), ")"
  )
)

C <- canCorPairs(form_check, vp_sampleinfo)

# Plot collinearity
png("figures/variance_partition_collinearity.png", units = "in", width = 8, height = 8, res = 300)
plotCorrMatrix(C)
dev.off()

# ─── Prepare count matrix ─────────────────────────────────────────────────────
# Keep only samples present in vp_sampleinfo
vp_samples <- rownames(vp_sampleinfo)
vp_counts  <- count_mat[, colnames(count_mat) %in% vp_samples, drop = FALSE]

# Reorder columns to match sampleinfo row order
vp_counts <- vp_counts[, vp_samples[vp_samples %in% colnames(vp_counts)], drop = FALSE]

# Sync sampleinfo to samples actually present in count matrix
vp_sampleinfo <- vp_sampleinfo[colnames(vp_counts), ]

cat("Samples in variance partition model:", ncol(vp_counts), "\n")
cat("Genes in variance partition model:  ", nrow(vp_counts), "\n")

# Normalize counts with edgeR CPM
dge        <- DGEList(counts = vp_counts)
dge        <- calcNormFactors(dge)
vp_cpm     <- cpm(dge, log = TRUE)  # log-CPM is standard for variancePartition

# Keep only genes with log-CPM > 1 in at least 10% of samples
keep <- rowSums(vp_cpm > 1) >= (0.1 * ncol(vp_cpm))
cat("Genes after filtering:", sum(keep), "\n")
vp_cpm_filt <- vp_cpm[keep, ]

# ─── Fit variance partition model (slow step) ────────────────────────────────
message("Fitting variance partition model — this may take a while...")

varPart <- fitExtractVarPartModel(vp_cpm_filt, form_check, vp_sampleinfo)

# Sort genes by median variance explained
vp_sorted <- sortCols(varPart)

# ─── Plot results ─────────────────────────────────────────────────────────────
# Violin plot per variable
vp_long <- pivot_longer(
  as.data.frame(vp_sorted),
  cols      = everything(),
  names_to  = "Variable",
  values_to = "VarianceExplained"
)
vp_long$Variable <- factor(vp_long$Variable, levels = colnames(vp_sorted))

# Clean up variable names for display
vp_long$Variable <- recode(vp_long$Variable,
                           "genes_contributing_to_80._of_reads" = "NG80",
                           "percentage_of_spliced_reads"        = "FSR",
                           "exonic_reads_minus_spike_ins"       = "FER",
                           "dataset_batch.y"                    = "Dataset",
                           "status"                             = "Phenotype",
                           "Residuals"                          = "Residuals",
                           "protein_coding_pct"                 = "Protein coding (%)",
                           "platelet"                           = "Platelet (%)"
)

p_vp <- ggplot(vp_long, aes(x = Variable, y = VarianceExplained, fill = Variable)) +
  geom_violin(scale = "width", alpha = 0.8) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = NULL,
    y = "Fraction of variance explained"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    text               = element_text(family = "Arial"),
    axis.title         = element_text(face = "bold", size = 12),
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y        = element_text(size = 10),
    legend.position    = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.8),
    plot.background    = element_rect(fill = "white", colour = "white")
  )

ggsave("figures/variance_partition_violin_no_dataset.png", p_vp,
       width = 8, height = 5, dpi = 600, device = ragg::agg_png)
ggsave("figures/variance_partition_violin_no_dataset.svg", p_vp,
       width = 8, height = 5, device = "svg")

# Save the variance partition table
write.table(vp_sorted, "tables/variance_partition_results.tsv",
            quote = FALSE, sep = "\t")

######################################################
# Variance Partition Analysis — Per Dataset
######################################################

library(BiocParallel)
library(parallel)
param <- SnowParam(detectCores() - 1, "SOCK", progressbar = TRUE)
register(param)

# Datasets to skip — parent-level collapsed groups that overlap with sub-batches
# to avoid double-counting samples
SKIP_DATASETS <- c("giraldez", "reggiardo", "block", "moufarrej", "roskams", "ibarra_plasma_cancer", "ibarra_plasma_non_cancer")
MIN_SAMPLES   <- 15  # minimum samples to attempt fitting

datasets <- levels(table_filtered$dataset_batch.y)
datasets <- datasets[!datasets %in% SKIP_DATASETS]

vp_results <- list()  # store varPart objects per dataset for later use

for (ds in datasets) {
  
  message("\n─── Processing dataset: ", ds, " ───")
  
  # ─── Subset to this dataset ──────────────────────────────────────────────
  ds_data <- filtered_df %>%
    filter(dataset_batch.y == ds) %>%
    select(sample_name, all_of(VP_NUMERIC), status) %>%
    filter(!is.na(genes_contributing_to_80._of_reads) &
             !is.na(percentage_of_spliced_reads) &
             !is.na(status)) %>%
    as.data.frame()
  
  row.names(ds_data) <- ds_data$sample_name
  
  # ─── Skip if too few samples ─────────────────────────────────────────────
  if (nrow(ds_data) < MIN_SAMPLES) {
    message("  Skipping — only ", nrow(ds_data), " samples (min: ", MIN_SAMPLES, ")")
    next
  }
  
  # ─── Check if status has more than 1 level ───────────────────────────────
  n_status <- length(unique(ds_data$status[!is.na(ds_data$status)]))
  has_status <- n_status > 1
  
  if (!has_status) {
    message("  Note: 'status' has only 1 level in this dataset — fitting without phenotype")
  }
  
  # ─── Scale numerics ──────────────────────────────────────────────────────
  ds_data$status <- as.factor(ds_data$status)
  for (num in VP_NUMERIC) {
    ds_data[[num]] <- scale(ds_data[[num]])
  }
  
  # ─── Build formula ───────────────────────────────────────────────────────
  if (has_status) {
    ds_form <- as.formula(
      paste0("~ ", paste(VP_NUMERIC, collapse = " + "), " + (1 | status)")
    )
  } else {
    ds_form <- as.formula(
      paste0("~ ", paste(VP_NUMERIC, collapse = " + "))
    )
  }
  
  # ─── Subset and normalize counts ─────────────────────────────────────────
  ds_samples <- rownames(ds_data)
  ds_counts  <- count_mat[, colnames(count_mat) %in% ds_samples, drop = FALSE]
  ds_counts  <- ds_counts[, ds_samples[ds_samples %in% colnames(ds_counts)], drop = FALSE]
  ds_data    <- ds_data[colnames(ds_counts), ]
  
  cat("  Samples:", ncol(ds_counts), "| ")
  
  dge_ds    <- DGEList(counts = ds_counts)
  dge_ds    <- calcNormFactors(dge_ds)
  cpm_ds    <- cpm(dge_ds, log = TRUE)
  
  # Filter lowly expressed genes
  keep      <- rowSums(cpm_ds > 1) >= (0.1 * ncol(cpm_ds))
  cpm_ds    <- cpm_ds[keep, ]
  cat("Genes after filtering:", nrow(cpm_ds), "\n")
  
  # ─── Fit model ───────────────────────────────────────────────────────────
  tryCatch({
    vp <- fitExtractVarPartModel(cpm_ds, ds_form, ds_data)
    vp_results[[ds]] <- vp
    
    # ─── Plot ──────────────────────────────────────────────────────────────
    vp_sorted_ds <- sortCols(vp)
    
    vp_long_ds <- pivot_longer(
      as.data.frame(vp_sorted_ds),
      cols      = everything(),
      names_to  = "Variable",
      values_to = "VarianceExplained"
    )
    vp_long_ds$Variable <- factor(vp_long_ds$Variable, levels = colnames(vp_sorted_ds))
    
    vp_long_ds$Variable <- recode(vp_long_ds$Variable,
                                  "genes_contributing_to_80._of_reads" = "NG80",
                                  "percentage_of_spliced_reads"        = "FSR",
                                  "exonic_reads_minus_spike_ins"       = "FER",
                                  "status"                             = "Phenotype",
                                  "Residuals"                          = "Residuals"
    )
    
    # Use dataset human-readable label as plot title if available
    ds_label <- if (!is.null(datasetsLabels[ds]) && !is.na(datasetsLabels[ds])) {
      datasetsLabels[ds]
    } else { ds }
    
    p_vp_ds <- ggplot(vp_long_ds, aes(x = Variable, y = VarianceExplained, fill = Variable)) +
      geom_violin(scale = "width", alpha = 0.8) +
      geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
      scale_fill_brewer(palette = "Set2") +
      scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
      labs(
        title = ds_label,
        x     = NULL,
        y     = "Fraction of variance explained"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        text               = element_text(family = "Arial"),
        axis.title         = element_text(face = "bold", size = 12),
        axis.text.x        = element_text(angle = 45, hjust = 1, size = 11),
        axis.text.y        = element_text(size = 10),
        plot.title         = element_text(face = "bold", size = 13),
        legend.position    = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.8),
        plot.background    = element_rect(fill = "white", colour = "white")
      )
    
    out_base <- paste0("figures/variance_partition_", ds)
    ggsave(paste0(out_base, ".png"), p_vp_ds,
           width = 6, height = 4, dpi = 600, device = ragg::agg_png)
    ggsave(paste0(out_base, ".svg"), p_vp_ds,
           width = 6, height = 4, device = "svg")
    
    message("  Saved: ", out_base)
    
  }, error = function(e) {
    message("  ERROR fitting model for ", ds, ": ", e$message)
  })
}

# ─── Optional: save all results ──────────────────────────────────────────────
saveRDS(vp_results, "tables/variance_partition_per_dataset.rds")
