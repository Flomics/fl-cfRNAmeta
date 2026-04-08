library(ggplot2)
library(dplyr)
library(ggpubr)
library(scales)
library(ggside)
library(ggnewscale)
library(grid)
library(jsonlite)
library(colorspace)
library(showtext)
library(svglite)
library("extrafont")
loadfonts()

# Read column names from text file
setwd("~/fl-cfRNAmeta/")
column_names <- c("read_number",
                  #"avg_input_read_length",
                  #"percentage_of_uniquely_mapped_reads",
                  "avg_mapped_read_length",
                  "mapped_percentage",
                  "exonic_percentage",                  
                  #"intronic_percentage",
                  "percentage_of_spliced_reads",
                  #"X.known_splice_junctions",           
                  #"read_coverage_uniformity_score",
                  #"junction_saturation_slope",
                  "median_insert_size",
                  "genes_contributing_to_80._of_reads",
                  "reads_mapping_sense_percentage",
                  "exonic_reads_minus_spike_ins",
                  "number_of_multimapped_reads",
                  "total_reads",
                  "number_of_uniquely_mapped_reads",
                  "spike_in_pct")
#"mt_rna_pct", "mt_rrna_pct", "mt_trna_pct", "misc_rna_pct", "protein_coding_pct", "lncrna_pct", "snrna_pct", "snorna_pct", "spike_in_pct", "other_rna_biotypes_pct")

data <- read.delim("tables/sampleinfo_all-batches.tsv", header = TRUE, sep = "\t", fileEncoding = "UTF-8")

# Load metadata
metadata <- read.delim("tables/cfRNA-meta_per_sample_metadata.tsv", header = TRUE, sep = "\t", fill = TRUE)


# Step 1: Filter merged_df to keep only samples that appear in metadata
metadata_subset <- metadata[, c("run", "dataset_batch")]

filtered_df <- data[data$sample_id %in% metadata$run, ]
filtered_df <- merge(filtered_df, metadata_subset, by.x = "sample_id", by.y = "run", all.x = TRUE)


removed_samples <- data[!(data$sample_id %in% metadata$run), ]
print(removed_samples$sample_id)


cat("Original merged_df rows:", nrow(data), "\n") #should be 2458 
cat("Filtered to samples in metadata:", nrow(filtered_df), "\n") # should be 2356, if it's not, make sure you have NOT removed all Flomics_1 samples due to a mismatch between the metadata names and the sampleinfo from snakeda names :)



biotype_data <- filtered_df[, 191:272] %>%
  select(-contains('_fc')) %>%
  mutate(across(everything(), ~ as.numeric(as.character(.))))


biotype_data$total <- rowSums(biotype_data)
biotype_data$percent_of_reads_mapping_to_spike_ins <- biotype_data$spike_in / biotype_data$total
biotype_data$percent_of_reads_mapping_to_spike_ins <- biotype_data$percent_of_reads_mapping_to_spike_ins * 100

exonic <- as.numeric(as.character(filtered_df$exonic))
spike_in <- as.numeric(as.character(filtered_df$spike_in))

filtered_df$exonic_reads_minus_spike_ins <- ifelse(
  is.na(exonic) & is.na(spike_in), NA,
  ifelse(is.na(exonic), 0, exonic) - ifelse(is.na(spike_in), 0, spike_in)
)


filtered_df$exonic_reads_minus_spike_ins <- (filtered_df$exonic_reads_minus_spike_ins / filtered_df$mapped_fragments)  * 100


# Keep specified columns
selected_columns <- NULL
selected_columns <- c( "sample_id", "sample_name", "sequencing_batch", "status", "dataset_batch.y", column_names)
filtered_data <- filtered_df[ ,selected_columns]

filtered_data$percent_of_reads_mapping_to_spike_ins <- biotype_data$percent_of_reads_mapping_to_spike_ins

result <- filtered_data %>%
  filter(percent_of_reads_mapping_to_spike_ins > 5) %>% 
  group_by(sequencing_batch) %>%     
  summarise(count = n())    

print(result)
# table_filtered <- filtered_data %>%
#   filter(percent_of_reads_mapping_to_spike_ins <= 5)
table_filtered <- filtered_data # NOT removing high spike-in samples

#table_filtered$log_genes_80 <- log(table_filtered$genes_contributing_to_80._of_reads)


num_datasets <- length(unique(table_filtered$dataset_batch.y))

mappings <- fromJSON("src/dataset_mappings.json")

datasetsLabels <- unlist(mappings$datasetsLabels)
core_order <- unlist(mappings$datasetVisualOrder)
datasetsPalette <- unlist(mappings$datasetsPalette)

table_filtered$dataset_batch.y <- factor(table_filtered$dataset_batch.y, levels = core_order)

column_names <- c(column_names, "percent_of_reads_mapping_to_spike_ins", "log_genes_80")

darken_color <- function(color, factor = 1.3) {
  rgb_col <- col2rgb(color) / 255
  darker_rgb <- pmin(rgb_col * (1 / factor), 1)
  rgb(darker_rgb[1], darker_rgb[2], darker_rgb[3])
}

datasetsOutlinePalette <- sapply(datasetsPalette, darken_color)

clean_label <- function(label) {
  label <- gsub("^X\\.", "", label)             
  label <- gsub("_", " ", label)                 
  label <- gsub("\\s+", " ", label)              
  label <- trimws(label)
  label <- tools::toTitleCase(label)              
  return(label)
}


bracket_df <- data.frame(
  xmin = c("block_150bp", "giraldez_phospho-rna-seq", "ibarra_buffy_coat", "reggiardo_bioivt", "moufarrej_site_1", "roskams_pilot"),
  xmax = c("block_300bp", "giraldez_standard", "ibarra_serum", "reggiardo_dls", "moufarrej_site_2", "roskams_validation"),
  label = c("Block", "Giráldez", "Ibarra", "Reggiardo", "Moufarrej", "Roskams-Hieter")
)

table_filtered$percent_of_multimapped_reads <- (table_filtered$number_of_multimapped_reads /table_filtered$total_reads)*100
column_names <- c(column_names, "percent_of_multimapped_reads")

table_filtered$percent_of_multimapped_reads_total_reads_mapped <- (table_filtered$number_of_multimapped_reads / (table_filtered$number_of_uniquely_mapped_reads + table_filtered$number_of_multimapped_reads ))*100
column_names <- c(column_names, "percent_of_multimapped_reads_total_reads_mapped")

#write.table(table_filtered, file="Qc_table_filtered.tsv", row.names = FALSE)
table_filtered$spike_in_pct <- table_filtered$spike_in_pct * 100

add_bottom_brackets <- function(p, bracket_df, factor_levels, y_base = -0.03, height = 0.015, col="black", lwd=0.8) {
  for (i in seq_len(nrow(bracket_df))) {
    x1 <- which(factor_levels == bracket_df$xmin[i])
    x2 <- which(factor_levels == bracket_df$xmax[i])
    if (length(x1) == 0 || length(x2) == 0) next
    
    offset <- 0.004
    x_start <- ((x1 - 1) / length(factor_levels)) + offset
    x_end   <- (x2 / length(factor_levels)) - offset
    
    bracket <- linesGrob(
      x = unit.c(unit(x_start, "npc"), unit(x_end, "npc")),
      y = unit(c(y_base, y_base), "npc"),
      gp = gpar(col = col, lwd = lwd)
    )
    
    verticals <- gList(
      linesGrob(
        x = unit.c(unit(x_start, "npc"), unit(x_start, "npc")),
        y = unit(c(y_base, y_base - height), "npc"),
        gp = gpar(col = col, lwd = lwd)
      ),
      linesGrob(
        x = unit.c(unit(x_end, "npc"), unit(x_end, "npc")),
        y = unit(c(y_base, y_base - height), "npc"),
        gp = gpar(col = col, lwd = lwd)
      )
    )
    
    p <- p + annotation_custom(grobTree(bracket, verticals))
  }
  return(p)
}


bpc_colors <- c(
  "cfDNA" =  "#AECAD9",
  "Custom" = "#D9BBAE",
  "Exome-based (EB)" = "#BBE0BB",
  "Whole RNA-Seq (oligo-dT pr.) (WRO)" = "#D0AED9",
  "Whole RNA-Seq (random pr.) (WRR)" = "#D9D6AE"
)

bpc_labels <- c(
  "Custom" = "Custom",
  "Exome-based (EB)" = "EB",
  "Whole RNA-Seq (oligo-dT pr.) (WRO)" = "WRO",
  "Whole RNA-Seq (random pr.) (WRR)" = "WRR",
  "cfDNA" = "cfDNA"
)

bpc_order <- c("Custom", "Exome-based (EB)", "Whole RNA-Seq (oligo-dT pr.) (WRO)", "Whole RNA-Seq (random pr.) (WRR)", "cfDNA")

ng_no_spike_ins <- read.delim("tables/genes_contributing_to_percentage_reads_no_spike_ins.tsv")

ng_no_spike_ins_table <- table_filtered %>%
  left_join(ng_no_spike_ins, by = c("sample_name" = "Sample"))

library(data.table)
raw_counts <- fread("tables/gene_raw_counts.tsv", sep = "\t")

# Get the sample columns that exist in both objects
sample_cols <- intersect(names(raw_counts), table_filtered$sample_name)
cat("Matched samples:", length(sample_cols), "\n")

# Subset the matrix to only matched samples
count_mat <- as.matrix(raw_counts[, ..sample_cols])  # data.table syntax for column selection

# Shannon Entropy per sample
shannon_entropy <- function(counts) {
  counts <- counts[counts > 0]          # remove zeros
  p <- counts / sum(counts)             # relative frequencies
  -sum(p * log2(p))                     # log base 2 (bits); use log() for nats
}

gini_index <- function(counts) {
  counts <- as.numeric(counts)
  counts <- counts[!is.na(counts) & counts >= 0]  # remove NAs before anything else
  n <- length(counts)
  if (n == 0 || sum(counts) == 0) return(NA)
  counts <- sort(counts)
  (2 * sum(seq_len(n) * counts)) / (n * sum(counts)) - (n + 1) / n
}


# Apply to each sample column
diversity_df <- data.frame(
  sample_name = sample_cols,
  shannon    = apply(count_mat, 2, shannon_entropy),
  gini       = apply(count_mat, 2, gini_index)
)

# Merge back into filtered_data
table_filtered <- merge(table_filtered, diversity_df, by = "sample_name", all.x = TRUE)

library(ggpubr)  # for stat_cor

# Ensure dataset factor levels and palette are set (as done earlier in your script)
table_filtered$dataset_batch.y <- factor(table_filtered$dataset_batch.y, levels = core_order)

adjusted_palette <- datasetsPalette[core_order]
names(adjusted_palette) <- core_order          # must match factor levels, i.e. raw IDs

# table_filtered$dataset_batch.y <- factor(
#   table_filtered$dataset_batch.y,
#   levels = core_order,
#   labels = datasetsLabels[core_order]
# )

# ─── Helper to build each scatterplot ────────────────────────────────────────
make_scatter <- function(data, x_var, y_var, x_label, y_label,
                         log_x = FALSE, log_y = FALSE, show_trend = TRUE, show_cor = TRUE) {
  
  p <- ggplot(data,
              aes(x = .data[[x_var]],
                  y = .data[[y_var]],
                  color = dataset_batch.y)) +
    geom_point(size = 2, alpha = 0.8) +
    scale_color_manual(
      values = adjusted_palette,
      labels = datasetsLabels[core_order],
      drop   = TRUE
    ) +
    labs(x = x_label, y = y_label, color = "Dataset") +
    theme_minimal(base_size = 13) +
    theme(
      text             = element_text(family = "Arial"),
      strip.text       = element_text(face = "bold"),
      axis.title       = element_text(face = "bold", size = 12),
      axis.text        = element_text(size = 10),
      legend.text      = element_text(size = 9),
      legend.title     = element_text(face = "bold"),
      legend.key.height = unit(0.5, "lines"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.8),
      panel.grid.minor.y = element_blank(),
      plot.background  = element_rect(fill = "white", colour = "white")
    ) +
    guides(color = guide_legend(ncol = 1))
  
  if (show_trend) p <- p + geom_smooth(method = "lm", se = FALSE, linewidth = 0.6,
                                       color = "grey50", alpha = 0.6)
  if (show_cor)   p <- p + stat_cor(method = "pearson", label.x.npc = "left",
                                    label.y.npc = "top", size = 3,
                                    aes(label = after_stat(r.label)), color = "black")
  if (log_x) p <- p + scale_x_continuous(trans = log10_trans())
  if (log_y) p <- p + scale_y_continuous(trans = log10_trans())
  
  return(p)
}



library(vegan)   # ecology diversity package — very well validated
library(ineq)    # inequality measures package

# ─── Package-based alternatives ──────────────────────────────────────────────

# Shannon via vegan::diversity() — standard in ecology/genomics literature
shannon_vegan <- function(counts) {
  counts <- as.numeric(counts)
  counts <- counts[!is.na(counts) & counts >= 0]
  vegan::diversity(counts, index = "shannon")  # uses natural log (nats), not log2
}

# Shannon via vegan using log2 
shannon_vegan_log2 <- function(counts) {
  counts <- as.numeric(counts)
  counts <- counts[!is.na(counts) & counts >= 0]
  vegan::diversity(counts, index = "shannon", base = 2) 
}

# Gini via ineq::Gini()
gini_ineq <- function(counts) {
  counts <- as.numeric(counts)
  counts <- counts[!is.na(counts) & counts >= 0]
  if (length(counts) == 0 || sum(counts) == 0) return(NA)
  ineq::Gini(counts)
}

# ─── Rebuild diversity_df with all implementations ───────────────────────────
diversity_df <- data.frame(
  sample_name       = sample_cols,
  shannon_manual    = apply(count_mat, 2, shannon_entropy),
  gini_manual       = apply(count_mat, 2, gini_index),
  # Package-based implementations
  shannon_vegan_nat = apply(count_mat, 2, shannon_vegan),       # nats
  shannon_vegan_log2= apply(count_mat, 2, shannon_vegan_log2),  # bits 
  gini_ineq         = apply(count_mat, 2, gini_ineq)
)

# Quick sanity check 
cat("Shannon manual vs vegan (log2) correlation: ",
    cor(diversity_df$shannon_manual, diversity_df$shannon_vegan_log2, use = "complete.obs"), "\n")
cat("Gini manual vs ineq correlation: ",
    cor(diversity_df$gini_manual, diversity_df$gini_ineq, use = "complete.obs"), "\n")

table_filtered <- merge(table_filtered, diversity_df, by = "sample_name", all.x = TRUE)

# ─── Plot 1: NG80 vs Shannon ──────────────────────────────────────────────────
p_ng80_shannon_vegan_nat <- make_scatter(
  data    = table_filtered,
  x_var   = "shannon_vegan_nat",
  y_var   = "genes_contributing_to_80._of_reads",
  x_label = "Shannon entropy (natural logarithm)",
  y_label = "NG80",
  log_y   = TRUE
)

ggsave("figures/ng80_vs_shannon_vegan_nat.png", p_ng80_shannon,
       width = 10, height = 6, dpi = 600, device = ragg::agg_png)
ggsave("figures/ng80_vs_shannon_vegan_nat.svg", p_ng80_shannon,
       width = 10, height = 6, device = "svg")

# ─── Plot 2: NG80 vs Gini ─────────────────────────────────────────────────────
p_ng80_gini <- make_scatter(
  data    = table_filtered,
  x_var   = "gini_ineq",
  y_var   = "genes_contributing_to_80._of_reads",
  x_label = "Gini index",
  y_label = "NG80",
  log_y   = TRUE,
  show_trend = FALSE,
  show_cor = FALSE
)

ggsave("figures/ng80_vs_gini_ineq.png", p_ng80_gini,
       width = 10, height = 6, dpi = 600, device = ragg::agg_png)
ggsave("figures/ng80_vs_gini_ineq.svg", p_ng80_gini,
       width = 10, height = 6, device = "svg")

# ─── Plot 3: Shannon vs Gini ──────────────────────────────────────────────────
p_shannon_gini <- make_scatter(
  data    = table_filtered,
  x_var   = "shannon_vegan_nat",
  y_var   = "gini_ineq",
  x_label = "Shannon entropy (natural logarithm)",
  y_label = "Gini index",
  show_trend = FALSE,
  show_cor = FALSE
)

ggsave("figures/shannon_vegan_nat_vs_gini_ineq.png", p_shannon_gini,
       width = 10, height = 6, dpi = 600, device = ragg::agg_png)
ggsave("figures/shannon_vegan_nat_vs_gini_ineq.svg", p_shannon_gini,
       width = 10, height = 6, device = "svg")

