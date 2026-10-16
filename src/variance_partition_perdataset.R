#!/usr/bin/env Rscript
# Variance partition analysis run separately within each study (dataset_short_name).
#
# Usage:
#   Rscript src/variance_partition_perdataset.R <sampleinfo> <metadata> <old_sampleinfo> \
#       <ng80_table> <tpm_matrix> <outdir> [min_samples] [cores]
#
#   sampleinfo      tables/sampleinfo_all-batches.tsv
#   metadata        tables/cfRNA-meta_per_sample_metadata.tsv
#   old_sampleinfo  tables/sampleinfo_external_and_internal_datasets.tsv (platelet fraction)
#   ng80_table      tables/genes_contributing_to_percentage_reads_no_spike_ins.tsv
#   tpm_matrix      tables/gene_tpm_norm.tsv
#   outdir          output directory (created if missing)
#   min_samples     studies with fewer samples are skipped (default 15)
#   cores           worker processes for the model fits (default: all but one)
#
# Writes per study: variance_partition_<study>.{png,svg,tsv} and
# variance_partition_collinearity_<study>.png; and variance_partition_perdataset_medians.tsv.

suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(variancePartition)
  library(BiocParallel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: Rscript src/variance_partition_perdataset.R <sampleinfo> <metadata> ",
       "<old_sampleinfo> <ng80_table> <tpm_matrix> <outdir> [min_samples] [cores]")
}
sampleinfo_file <- args[1]
metadata_file   <- args[2]
old_si_file     <- args[3]
ng80_file       <- args[4]
tpm_file        <- args[5]
outdir          <- args[6]
min_samples     <- if (length(args) >= 7) as.integer(args[7]) else 15
cores           <- if (length(args) >= 8) as.integer(args[8]) else max(1, parallel::detectCores() - 1)

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
param <- MulticoreParam(cores)

# ─── Sample table ────────────────────────────────────────────────────────────
si <- fread(sampleinfo_file, sep = "\t", quote = "", data.table = FALSE)
md <- fread(metadata_file, sep = "\t", data.table = FALSE)

# The sampleinfo has samples the matrix does not; the metadata is the analysed set.
df <- md %>%
  select(run, dataset_short_name, dataset_batch, broad_protocol_category) %>%
  inner_join(si %>% select(sample_id, sample_name, phenotype, percentage_of_spliced_reads,
                           mapped_fragments, qualimap_exonic, qualimap_intronic,
                           qualimap_intergenic, spike_in),
             by = c("run" = "sample_id"))

# FER as in boxplots_fig2.R
exonic   <- as.numeric(df$qualimap_exonic)
spike_in <- as.numeric(df$spike_in)
df$exonic_reads_minus_spike_ins <- ifelse(
  is.na(exonic) & is.na(spike_in), NA,
  ifelse(is.na(exonic), 0, exonic) - ifelse(is.na(spike_in), 0, spike_in)
) / (df$qualimap_exonic + df$qualimap_intronic + df$qualimap_intergenic) * 100

ng80 <- fread(ng80_file, sep = "\t", data.table = FALSE)
df$ng80 <- ng80[["number_of_genes_contributing_to_80%_of_reads"]][match(df$sample_name, ng80$Sample)]

# Platelet fractions only exist in the old run's sampleinfo, which uses old sample
# names: link on run, then sample_id, then current sample_name on old sample_id.
old <- fread(old_si_file, sep = "\t", quote = "", data.table = FALSE)
platelet_by_run <- setNames(old$platelet, old$run)
platelet_by_sid <- setNames(old$platelet, old$sample_id)
df$platelet <- coalesce(platelet_by_run[df$run],
                        platelet_by_sid[df$run],
                        platelet_by_sid[df$sample_name])
cat("Samples with a platelet fraction:", sum(!is.na(df$platelet)), "of", nrow(df), "\n")

# Same collapsing as variance_partition.R; later lines overwrite earlier ones.
simplify_phenotype <- function(phenotype) {
  out <- rep(NA_character_, length(phenotype))
  out[phenotype == "healthy"] <- "healthy"
  out[is.na(phenotype) | phenotype == ""] <- "missing"
  rules <- c(
    "Acute Myeloid Leukemia"                       = "cancer",
    "Alzheimers disease"                           = "non-cancer disease",
    "Chronic hepatitis B"                          = "non-cancer disease",
    "Chronic kidney failure EPO-treated"           = "non-cancer disease",
    "Cirrhosis"                                    = "non-cancer disease",
    "Colorectal cancer"                            = "cancer",
    "Diffuse large B-cell lymphoma"                = "cancer",
    "Diverticulitis"                               = "non-cancer disease",
    "Esophagus cancer"                             = "cancer",
    "G-CSF-treated healthy donors"                 = "healthy",
    "Healthy"                                      = "healthy",
    "Healthy pregnant woman"                       = "healthy",
    "Healthy pregnant woman who delivered preterm" = "non-cancer disease",
    "Liver cancer"                                 = "cancer",
    "Lung cancer"                                  = "cancer",
    "Multiple myeloma"                             = "cancer",
    "Nonalcoholic fatty liver disease"             = "non-cancer disease",
    "Nonalcoholic steatohepatitis"                 = "non-cancer disease",
    "Pancreatic cancer"                            = "cancer",
    "Pre-cancerous condition: cirrhosis"           = "non-cancer disease",
    "Pre-cancerous condition: MGUS"                = "non-cancer disease",
    "Pre-eclampsia"                                = "non-cancer disease",
    "Primary mediastinal B-cell lymphoma"          = "cancer",
    "Stomach cancer"                               = "cancer"
  )
  for (pattern in names(rules)) out[grep(pattern, phenotype)] <- rules[[pattern]]
  out
}
df$simple_phenotype <- simplify_phenotype(df$phenotype)

VP_NUMERIC     <- c("ng80", "percentage_of_spliced_reads", "exonic_reads_minus_spike_ins",
                    "platelet", "mapped_fragments")
VP_CATEGORICAL <- c("dataset_batch", "simple_phenotype", "broad_protocol_category")

VAR_LABELS <- c(
  ng80                         = "NG80",
  percentage_of_spliced_reads  = "FSR",
  exonic_reads_minus_spike_ins = "FER",
  platelet                     = "Platelet (%)",
  mapped_fragments             = "Mapped fragments",
  dataset_batch                = "Batch",
  simple_phenotype             = "Phenotype",
  broad_protocol_category      = "Broad Protocol Category (BPC)",
  Residuals                    = "Residuals"
)
var_cols <- setNames(colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(VAR_LABELS)),
                     VAR_LABELS)

source("src/load_dataset_mappings.R")
ds_labels <- load_dataset_mappings("src/dataset_mappings.json")$datasetsLabels
study_label <- function(study) {
  if (!is.na(ds_labels[study])) return(unname(ds_labels[study]))
  tools::toTitleCase(study)
}

# ─── TPM matrix ──────────────────────────────────────────────────────────────
tpm <- fread(tpm_file, sep = "\t", data.table = FALSE)
rownames(tpm) <- tpm$gene_id
tpm <- as.matrix(tpm[, intersect(colnames(tpm), df$sample_name)])
cat("Samples in both matrix and metadata:", ncol(tpm), "\n")

# ─── One model per study ─────────────────────────────────────────────────────
medians <- list()

for (study in sort(unique(df$dataset_short_name))) {
  message("\n─── ", study, " ───")
  d <- df %>% filter(dataset_short_name == study, sample_name %in% colnames(tpm))

  # A variable entirely missing in this study is dropped, not its samples.
  numeric_vars <- VP_NUMERIC[sapply(VP_NUMERIC, function(v) any(!is.na(d[[v]])))]
  n_before <- nrow(d)
  d <- d[complete.cases(d[, numeric_vars]), ]
  if (nrow(d) < n_before) message("  Dropped ", n_before - nrow(d), " samples with missing values")

  if (nrow(d) < min_samples) {
    message("  Skipping: ", nrow(d), " samples (min ", min_samples, ")")
    next
  }

  # Within one study most of these are constant, and a constant term cannot be fitted.
  numeric_vars <- numeric_vars[sapply(numeric_vars, function(v) sd(d[[v]]) > 0)]
  cat_vars <- VP_CATEGORICAL[sapply(VP_CATEGORICAL, function(v) length(unique(d[[v]])) > 1)]
  message("  Samples: ", nrow(d), " | numeric: ", paste(numeric_vars, collapse = ", "),
          " | categorical: ", if (length(cat_vars)) paste(cat_vars, collapse = ", ") else "none")

  d$ng80 <- log(d$ng80)
  for (v in numeric_vars) d[[v]] <- as.numeric(scale(d[[v]]))
  for (v in cat_vars) d[[v]] <- droplevels(as.factor(d[[v]]))
  rownames(d) <- d$sample_name

  terms <- c(numeric_vars, if (length(cat_vars)) paste0("(1 | ", cat_vars, ")"))
  form  <- as.formula(paste("~", paste(terms, collapse = " + ")))

  C <- canCorPairs(as.formula(paste("~", paste(c(numeric_vars, cat_vars), collapse = " + "))), d)
  png(file.path(outdir, paste0("variance_partition_collinearity_", study, ".png")),
      units = "in", width = 8, height = 8, res = 300)
  plotCorrMatrix(C)
  dev.off()

  # TPM is per-sample, so this study's columns are its own TPM; the expression
  # filter is recomputed on them.
  x <- tpm[, d$sample_name, drop = FALSE]
  x <- log2(x[rowSums(x > 1) >= 0.1 * ncol(x), , drop = FALSE] + 1)
  message("  Genes after filtering: ", nrow(x))

  vp <- tryCatch(
    sortCols(fitExtractVarPartModel(x, form, d, BPPARAM = param)),
    error = function(e) { message("  ERROR: ", conditionMessage(e)); NULL }
  )
  if (is.null(vp)) next

  vp_df <- as.data.frame(vp)
  write.table(cbind(gene_id = rownames(vp_df), vp_df),
              file.path(outdir, paste0("variance_partition_", study, ".tsv")),
              quote = FALSE, sep = "\t", row.names = FALSE)
  medians[[study]] <- data.frame(study = study, n_samples = nrow(d), n_genes = nrow(x),
                                 variable = colnames(vp_df),
                                 median_variance = apply(vp_df, 2, median))

  vp_long <- pivot_longer(vp_df, everything(), names_to = "Variable", values_to = "VarianceExplained")
  vp_long$Variable <- factor(VAR_LABELS[vp_long$Variable], levels = VAR_LABELS[colnames(vp_df)])

  p <- ggplot(vp_long, aes(x = Variable, y = VarianceExplained, fill = Variable)) +
    geom_violin(scale = "width", alpha = 0.8) +
    geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
    scale_fill_manual(values = var_cols) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(title = paste0(study_label(study), " (n = ", nrow(d), ")"),
         x = NULL, y = "Fraction of variance explained") +
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
      panel.grid.major.y = element_line(linewidth = 0.8, colour = "grey90"),
      plot.background    = element_rect(fill = "white", colour = "white")
    )

  out_base <- file.path(outdir, paste0("variance_partition_", study))
  ggsave(paste0(out_base, ".png"), p, width = 6, height = 4, dpi = 600, device = ragg::agg_png)
  ggsave(paste0(out_base, ".svg"), p, width = 6, height = 4, device = svglite::svglite)
}

write.table(bind_rows(medians), file.path(outdir, "variance_partition_perdataset_medians.tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE)
