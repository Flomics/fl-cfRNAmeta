######################################################
# Variance Partition Analysis
# Needs to run after the boxplots_fig2.R script
######################################################

suppressMessages(library("variancePartition"))
suppressMessages(library("edgeR"))

# ─── Variables to test ───────────────────────────────────────────────────────
VP_NUMERIC <- c(
  "genes_contributing_to_80._of_reads",   # NG80
  "percentage_of_spliced_reads",           # FSR
  "exonic_reads_minus_spike_ins"           # FER (optional, remove if not needed)
)

VP_CATEGORICAL <- c(
  "dataset_batch.y",      # dataset of origin — likely major driver
  "status"                # phenotype (case/control)
)

# ─── Prepare sampleinfo ──────────────────────────────────────────────────────
# Start from table_filtered which already has all the QC metrics
vp_sampleinfo <- table_filtered %>%
  select(sample_name, all_of(VP_NUMERIC), all_of(VP_CATEGORICAL)) %>%
  filter(!is.na(genes_contributing_to_80._of_reads) &
           !is.na(percentage_of_spliced_reads) &
           !is.na(status)) %>%
  as.data.frame()

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
                           "sequencing_batch"                   = "Seq. batch",
                           "status"                             = "Phenotype",
                           "Residuals"                          = "Residuals"
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

ggsave("figures/variance_partition_violin.png", p_vp,
       width = 8, height = 5, dpi = 600, device = ragg::agg_png)
ggsave("figures/variance_partition_violin.svg", p_vp,
       width = 8, height = 5, device = "svg")

# Save the variance partition table
write.table(vp_sorted, "tables/variance_partition_results.tsv",
            quote = FALSE, sep = "\t")