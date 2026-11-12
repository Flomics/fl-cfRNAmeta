#!/usr/bin/env Rscript
# Variance partition within each study (dataset_short_name).
#
# Writes per study a violin plot, the variance of each gene and a collinearity matrix, plus
# variance_partition_perdataset_{medians,summary}.tsv with all studies.
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
#   min_samples     default 15
#   cores           default: all but one

suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(variancePartition)
  library(BiocParallel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) stop("see the usage block at the top of the script")
outdir      <- args[6]
min_samples <- if (length(args) >= 7) as.integer(args[7]) else 15
cores       <- if (length(args) >= 8) as.integer(args[8]) else max(1, parallel::detectCores() - 1)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

si  <- fread(args[1], sep = "\t", quote = "", data.table = FALSE)
md  <- fread(args[2], sep = "\t", data.table = FALSE)
old <- fread(args[3], sep = "\t", quote = "", data.table = FALSE)
ng80 <- fread(args[4], sep = "\t", data.table = FALSE)

df <- md %>%
  select(run, dataset_short_name, dataset_batch, broad_protocol_category,
         collection_center_short_name) %>%
  inner_join(si, by = c("run" = "sample_id"), suffix = c("", ".si")) %>%
  mutate(
    # FER as in boxplots_fig2.R
    exonic_reads_minus_spike_ins = (coalesce(qualimap_exonic, 0) - coalesce(spike_in, 0)) /
      (qualimap_exonic + qualimap_intronic + qualimap_intergenic) * 100,
    ng80 = ng80[["number_of_genes_contributing_to_80%_of_reads"]][match(sample_name, ng80$Sample)],
    # the old sampleinfo has the old sample names
    platelet = coalesce(old$platelet[match(run, old$run)],
                        old$platelet[match(run, old$sample_id)],
                        old$platelet[match(sample_name, old$sample_id)])
  )

# within a study, the study's own labels rather than healthy / cancer / non-cancer
df$phenotype[df$phenotype == ""] <- NA
df$collection_centre <- gsub(",(\\S)", ", \\1", trimws(df$collection_center_short_name))
df$collection_centre[is.na(df$collection_centre) | df$collection_centre %in% c("", "Unknown")] <- "Unspecified"

numeric_all     <- c("ng80", "percentage_of_spliced_reads", "exonic_reads_minus_spike_ins",
                     "platelet", "mapped_fragments")
# order matters: of two terms that split the samples identically, the first one is kept
categorical_all <- c("phenotype", "collection_centre", "dataset_batch", "broad_protocol_category")

labels <- c(ng80 = "NG80", percentage_of_spliced_reads = "FSR",
            exonic_reads_minus_spike_ins = "FER", platelet = "Platelet (%)",
            mapped_fragments = "Mapped fragments", dataset_batch = "Batch",
            phenotype = "Phenotype", collection_centre = "Collection centre",
            broad_protocol_category = "Broad Protocol Category (BPC)",
            Residuals = "Residuals")
colours <- setNames(colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(labels)), labels)

source("src/load_dataset_mappings.R")
ds_labels <- load_dataset_mappings("src/dataset_mappings.json")$datasetsLabels

tpm <- fread(args[5], sep = "\t", data.table = FALSE)
rownames(tpm) <- tpm$gene_id
tpm <- as.matrix(tpm[, intersect(colnames(tpm), df$sample_name)])

medians <- list()

for (study in sort(unique(df$dataset_short_name))) {
  d <- df %>% filter(dataset_short_name == study, sample_name %in% colnames(tpm))
  numeric_vars <- numeric_all[sapply(numeric_all, function(v) any(!is.na(d[[v]])))]
  d <- d[complete.cases(d[, numeric_vars]) & !is.na(d$phenotype), ]
  if (nrow(d) < min_samples) {
    message(study, ": skipped, ", nrow(d), " samples")
    next
  }

  # terms that are constant within the study cannot be fitted
  numeric_vars <- numeric_vars[sapply(numeric_vars, function(v) sd(d[[v]]) > 0)]
  cat_vars <- categorical_all[sapply(categorical_all, function(v) n_distinct(d[[v]]) > 1)]

  # perfectly confounded terms cannot be told apart: keep one and record it
  same_split <- function(a, b) n_distinct(paste(a, b)) == n_distinct(a) && n_distinct(a) == n_distinct(b)
  kept <- character()
  confounded <- character()
  for (v in cat_vars) {
    hit <- Filter(function(k) same_split(d[[k]], d[[v]]), kept)
    if (length(hit)) confounded <- c(confounded, paste(v, "=", hit[1])) else kept <- c(kept, v)
  }
  cat_vars <- kept
  confounded <- paste(confounded, collapse = "; ")
  message(study, ": ", nrow(d), " samples, ", paste(c(numeric_vars, cat_vars), collapse = " + "),
          if (nzchar(confounded)) paste0(" (dropped, confounded: ", confounded, ")"))

  d$ng80 <- log(d$ng80)
  d[numeric_vars] <- lapply(d[numeric_vars], function(x) as.numeric(scale(x)))
  d[cat_vars] <- lapply(d[cat_vars], factor)
  rownames(d) <- d$sample_name

  png(file.path(outdir, paste0("variance_partition_collinearity_", study, ".png")),
      units = "in", width = 8, height = 8, res = 300)
  plotCorrMatrix(canCorPairs(reformulate(c(numeric_vars, cat_vars)), d))
  dev.off()

  # TPM is normalised per sample, so a study's columns are its own TPM
  x <- tpm[, d$sample_name]
  x <- log2(x[rowSums(x > 1) >= 0.1 * ncol(x), ] + 1)

  form <- reformulate(c(numeric_vars, if (length(cat_vars)) paste0("(1 | ", cat_vars, ")")))
  vp <- tryCatch(fitExtractVarPartModel(x, form, d, BPPARAM = MulticoreParam(cores)),
                 error = function(e) { message(study, ": ", conditionMessage(e)); NULL })
  if (is.null(vp)) next
  vp <- as.data.frame(sortCols(vp))

  write.table(cbind(gene_id = rownames(vp), vp),
              file.path(outdir, paste0("variance_partition_", study, ".tsv")),
              quote = FALSE, sep = "\t", row.names = FALSE)
  medians[[study]] <- data.frame(study = study, n_samples = nrow(d), n_genes = nrow(x),
                                 variable = colnames(vp), median_variance = apply(vp, 2, median),
                                 genes_over_10pct = colMeans(vp > 0.1),
                                 confounded = confounded)

  long <- pivot_longer(vp, everything(), names_to = "Variable", values_to = "VarianceExplained")
  long$Variable <- factor(labels[long$Variable], levels = labels[colnames(vp)])
  title <- if (is.na(ds_labels[study])) tools::toTitleCase(study) else ds_labels[[study]]

  p <- ggplot(long, aes(x = Variable, y = VarianceExplained, fill = Variable)) +
    geom_violin(scale = "width", alpha = 0.8) +
    geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
    scale_fill_manual(values = colours) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(title = paste0(title, " (n = ", nrow(d), ")"), x = NULL,
         y = "Fraction of variance explained") +
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

  out <- file.path(outdir, paste0("variance_partition_", study))
  ggsave(paste0(out, ".png"), p, width = 6, height = 4, dpi = 600, device = ragg::agg_png)
  ggsave(paste0(out, ".svg"), p, width = 6, height = 4, device = svglite::svglite)
}

medians <- bind_rows(medians)
write.table(medians, file.path(outdir, "variance_partition_perdataset_medians.tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE)

# one row per study, in percent: the median alone hides the genes a term does explain
summary <- medians %>%
  mutate(variable = factor(variable, levels = c(categorical_all, numeric_all, "Residuals")),
         median_pct = round(100 * median_variance, 1),
         genes_over_10pct = round(100 * genes_over_10pct, 1)) %>%
  arrange(variable) %>%
  pivot_wider(id_cols = c(study, n_samples, n_genes, confounded), names_from = variable,
              values_from = c(median_pct, genes_over_10pct),
              names_glue = "{variable}_{.value}", names_vary = "slowest") %>%
  relocate(confounded, .after = last_col()) %>%
  arrange(desc(n_samples))
write.table(summary, file.path(outdir, "variance_partition_perdataset_summary.tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE, na = "NA")
