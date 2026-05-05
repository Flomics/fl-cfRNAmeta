library(dplyr)
library(tidyr)
library(data.table)

######################################################
# Variance Partition Analysis
# Needs to run after the boxplots_fig2.R script
######################################################

suppressMessages(library("variancePartition"))
suppressMessages(library("edgeR"))

setwd("~/fl-cfRNAmeta")

platelet_info <- read.delim("tables/sampleinfo_external_and_internal_datasets.tsv", header = TRUE, sep = "\t", fileEncoding = "UTF-8")

# Add platelet info to filtered_df so it's available for vp_sampleinfo
filtered_df <- filtered_df %>%
  left_join(
    platelet_info %>% select(sample_name, platelet),
    by = "sample_name"
  ) %>%
  left_join(
    metadata %>% select(run, broad_protocol_category),
    by = c("sample_id" = "run")
    )

filtered_df$simple_phenotype <- NA
filtered_df$simple_phenotype[filtered_df$phenotype == "healthy"] <- "healthy"
filtered_df$simple_phenotype[is.na(filtered_df$phenotype)] <- "missing"
filtered_df$simple_phenotype[filtered_df$phenotype == ""] <- "missing"
filtered_df$simple_phenotype[grep("Acute Myeloid Leukemia",filtered_df$phenotype)] <- "cancer"
filtered_df$simple_phenotype[grep("Alzheimers disease",filtered_df$phenotype)] <- "non-cancer disease"
filtered_df$simple_phenotype[grep("Chronic hepatitis B",filtered_df$phenotype)] <- "non-cancer disease"
filtered_df$simple_phenotype[grep("Chronic kidney failure EPO-treated",filtered_df$phenotype)] <- "non-cancer disease"
filtered_df$simple_phenotype[grep("Cirrhosis",filtered_df$phenotype)] <- "non-cancer disease"
filtered_df$simple_phenotype[grep("Colorectal cancer",filtered_df$phenotype)] <- "cancer"
filtered_df$simple_phenotype[grep("Diffuse large B-cell lymphoma",filtered_df$phenotype)] <- "cancer"
filtered_df$simple_phenotype[grep("Diverticulitis",filtered_df$phenotype)] <- "non-cancer disease"
filtered_df$simple_phenotype[grep("Esophagus cancer",filtered_df$phenotype)] <- "cancer"
filtered_df$simple_phenotype[grep("G-CSF-treated healthy donors",filtered_df$phenotype)] <- "healthy"
filtered_df$simple_phenotype[grep("Healthy",filtered_df$phenotype)] <- "healthy"
filtered_df$simple_phenotype[grep("Healthy pregnant woman",filtered_df$phenotype)] <- "healthy"
filtered_df$simple_phenotype[grep("Healthy pregnant woman who delivered preterm",filtered_df$phenotype)] <- "non-cancer disease"
filtered_df$simple_phenotype[grep("Liver cancer",filtered_df$phenotype)] <- "cancer"
filtered_df$simple_phenotype[grep("Lung cancer",filtered_df$phenotype)] <- "cancer"
filtered_df$simple_phenotype[grep("Multiple myeloma",filtered_df$phenotype)] <- "cancer"
filtered_df$simple_phenotype[grep("Nonalcoholic fatty liver disease",filtered_df$phenotype)] <- "non-cancer disease"
filtered_df$simple_phenotype[grep("Nonalcoholic steatohepatitis",filtered_df$phenotype)] <- "non-cancer disease"
filtered_df$simple_phenotype[grep("Pancreatic cancer",filtered_df$phenotype)] <- "cancer"
filtered_df$simple_phenotype[grep("Pre-cancerous condition: cirrhosis",filtered_df$phenotype)] <- "non-cancer disease"
filtered_df$simple_phenotype[grep("Pre-cancerous condition: MGUS",filtered_df$phenotype)] <- "non-cancer disease"
filtered_df$simple_phenotype[grep("Pre-eclampsia",filtered_df$phenotype)] <- "non-cancer disease"
filtered_df$simple_phenotype[grep("Primary mediastinal B-cell lymphoma",filtered_df$phenotype)] <- "cancer"
filtered_df$simple_phenotype[grep("Stomach cancer",filtered_df$phenotype)] <- "cancer"



# # ─── Datasets to include in the analysis ─────────────────────────────────────
# DATASETS_TO_INCLUDE <- c(
#   "chen", "decruyenaere", "flomics_2", "moufarrej_site_1", "moufarrej_site_2", 
#   "roskams_pilot", "roskams_validation", "tao", "zhu"
# )
# 
# filtered_df <- filtered_df %>%
#   filter(dataset_batch.y %in% DATASETS_TO_INCLUDE)
# 
# cat("Datasets included:", paste(DATASETS_TO_INCLUDE, collapse = ", "), "\n")
# cat("Samples after dataset filter:", nrow(filtered_df), "\n")

# ─── Variables to test ───────────────────────────────────────────────────────
VP_NUMERIC <- c(
  "genes_contributing_to_80._of_reads",   # NG80
  "percentage_of_spliced_reads",           # FSR
  "protein_coding_pct",                     # biotype composition
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
  "platelet",
  "mapped_fragments",
  "read_number"
)

VP_CATEGORICAL <- c(
  "dataset_batch.y",      # dataset
  "simple_phenotype",           # phenotype
  "read_length",
  "plasma_tubes",
  "biomaterial",
  "nucleic_acid_type",
  "rna_extraction_kit_short_name",
  "dnase",
  "library_prep_kit_short_name",
  "library_selection",
  "cdna_library_type",
  "centrifugation_step_1",
  "centrifugation_step_2",
  "broad_protocol_category.y"
)


VP_NUMERIC <- c(
  "genes_contributing_to_80._of_reads",   # NG80
  "percentage_of_spliced_reads",           # FSR
  "exonic_reads_minus_spike_ins",          # FER
  "platelet.y",
  "mapped_fragments",
  "read_number"
)

VP_CATEGORICAL <- c(
  "dataset_batch.y",      # dataset
  "simple_phenotype",           # phenotype
  "broad_protocol_category"
)

# ─── Prepare sampleinfo ──────────────────────────────────────────────────────
# Start from table_filtered which already has all the QC metrics
# vp_sampleinfo <- table_filtered %>%
#   select(sample_name, all_of(VP_NUMERIC), all_of(VP_CATEGORICAL)) %>%
#   filter(!is.na(genes_contributing_to_80._of_reads) &
#            !is.na(percentage_of_spliced_reads) &
#            !is.na(status)) %>%
#   as.data.frame()

vp_sampleinfo <- filtered_df %>%
  select(sample_name, dataset_batch.y, all_of(VP_NUMERIC), all_of(VP_CATEGORICAL)) %>%
  # Only require the complete variables — leave NAs for borderline ones
  filter(
    !is.na(genes_contributing_to_80._of_reads),
    !is.na(percentage_of_spliced_reads),
    !is.na(simple_phenotype)
  ) %>%
  as.data.frame()

# Report how many samples remain
cat("Samples after NA filtering:", nrow(vp_sampleinfo), "\n")
cat("Samples dropped:", nrow(filtered_df) - nrow(vp_sampleinfo), "\n")

row.names(vp_sampleinfo) <- vp_sampleinfo$sample_name

# ─── Clean categorical variables ─────────────────────────────────────────────
# Sentinel strings that should be treated as missing
INVALID_LEVELS <- c("Unspecified", "unspecified", "None", "none", "NA", "N/A", "n/a", "")

for (cat in VP_CATEGORICAL) {
  x <- trimws(as.character(vp_sampleinfo[[cat]]))
  invalid <- x %in% INVALID_LEVELS | grepl("^\\s*$", x) | is.na(x)
  if (any(invalid)) {
    cat("  Recoding", sum(invalid), "invalid values in", cat, "-> 'Other'\n")
    x[invalid] <- "Other"
  }
  vp_sampleinfo[[cat]] <- x
}

# Make categorical variables factors
for (cat in VP_CATEGORICAL) {
  vp_sampleinfo[[cat]] <- droplevels(as.factor(vp_sampleinfo[[cat]]))
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

# Log-transform NG80 before scaling
vp_sampleinfo$genes_contributing_to_80._of_reads <- log(vp_sampleinfo$genes_contributing_to_80._of_reads)

# Scale numeric variables (important for variance partition)
for (num in VP_NUMERIC) {
  vp_sampleinfo[[num]] <- as.numeric(scale(vp_sampleinfo[[num]]))
}

# ─── Check collinearity between metadata variables ────────────────────────────
# Flat formula for canCorPairs (does not accept random effects)
form_canCor <- as.formula(
  paste0("~ ", paste(c(VP_NUMERIC, VP_CATEGORICAL), collapse = " + "))
)
C <- canCorPairs(form_canCor, vp_sampleinfo)

# Model formula for VPA: numerics as fixed effects, categoricals as random effects
form_check <- as.formula(
  paste0(
    "~ ",
    paste(VP_NUMERIC, collapse = " + "),
    " + (1 | ", paste(VP_CATEGORICAL, collapse = ") + (1 | "), ")"
  )
)

# Plot collinearity
png("figures/variance_partition_collinearity_figure1_filtered.png", units = "in", width = 8, height = 8, res = 300)
plotCorrMatrix(C)
dev.off()



# ─── Prepare count matrix ─────────────────────────────────────────────────────
# Keep only samples present in vp_sampleinfo
vp_samples <- rownames(vp_sampleinfo)
vp_counts  <- count_mat[, colnames(count_mat) %in% vp_samples, drop = FALSE]
tmp <- fread("tables/gene_tpm_norm.tsv")
sample_cols <- intersect(names(tmp), table_filtered$sample_name)
cat("Matched samples:", length(sample_cols), "\n")

vp_counts <- as.matrix(tmp[, ..sample_cols])
rm(tmp); gc()

# Subset and reorder to vp_samples
vp_counts <- vp_counts[, vp_samples[vp_samples %in% colnames(vp_counts)], drop = FALSE]

# Sync sampleinfo to samples actually present in matrix
vp_sampleinfo <- vp_sampleinfo[colnames(vp_counts), ]

cat("Samples in variance partition model:", ncol(vp_counts), "\n")
cat("Genes in variance partition model:  ", nrow(vp_counts), "\n")

# Filter: keep genes with TPM > 1 in at least 10% of samples, then log2-transform
keep        <- rowSums(vp_counts > 1) >= (0.1 * ncol(vp_counts))
cat("Genes after filtering:", sum(keep), "\n")
gc()
vp_tpm_filt <- log2(vp_counts[keep, ] + 1)
rm(vp_counts); gc()

# ─── Parallelisation (register before the slow model fit) ────────────────────
library(BiocParallel)
library(parallel)
param <- MulticoreParam(detectCores() - 1, progressbar = TRUE)
register(param)

# ─── Fit variance partition model (slow step) ────────────────────────────────
message("Fitting variance partition model — this may take a while...")
gc()
varPart <- fitExtractVarPartModel(vp_tpm_filt, form_check, vp_sampleinfo, BPPARAM = param)

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
                           "platelet.y"                           = "Platelet (%)",
                           "broad_protocol_category" = "Broad Protocol Category (BPC)",
                           "mapped_fragments" = "Mapped fragments",
                           "cdna_library_type" = "cDNA Library type",
                           "simple_phenotype" = "Phenotype",
                           "read_number" = "Read number",
                           "read_length" = "Read length",
                           "plasma_tubes" = "Plasma tubes",
                           "biomaterial" = "Biomaterial",
                           "nucleic_acid_type" = "Nucleic acid type",
                           "rna_extraction_kit_short_name" = "RNA extraction kit",
                           "dnase" = "DNase",
                           "library_prep_kit_short_name" = "Library prep kit",
                           "library_selection" = "Library selection",
                           "centrifugation_step_1" = "Centrifugation step 1",
                           "centrifugation_step_2"  = "Centrifugation step 2"
)

n_vars   <- length(levels(vp_long$Variable))
vp_cols  <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n_vars)

p_vp <- ggplot(vp_long, aes(x = Variable, y = VarianceExplained, fill = Variable)) +
  geom_violin(scale = "width", alpha = 0.8) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
  scale_fill_manual(values = vp_cols) +
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

ggsave("figures/variance_partition_violin_figure_1_filtered.png", p_vp,
       width = 8, height = 5, dpi = 600, device = ragg::agg_png)
 ggsave("figures/variance_partition_violin_logng80.svg", p_vp,
       width = 8, height = 5, device = "svg")

# Save the variance partition table
write.table(vp_sorted, "tables/variance_partition_results.tsv",
            quote = FALSE, sep = "\t")

######################################################
# Variance Partition Analysis — Per Dataset
######################################################

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
  ds_data$genes_contributing_to_80._of_reads <- log(ds_data$genes_contributing_to_80._of_reads)
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
