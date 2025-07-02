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
#font_add(family="Arial", regular = "/usr/share/fonts/truetype/msttcorefonts/Arial.ttf")
#showtext_opts(dpi = 600)  # MUST come before showtext_auto()
#showtext_auto()
#theme_set(theme_classic(base_family = "Arial"))
library("extrafont")


# Read column names from text file
setwd("~/fl-cfRNAmeta/")
column_names <- c("read_number",
                  "avg_input_read_length",
                  "percentage_of_uniquely_mapped_reads",
                  "avg_mapped_read_length",
                  "mapped_percentage",
                  "exonic_percentage",                  
                  "intronic_percentage",
                  "percentage_of_spliced_reads",
                  "X.known_splice_junctions",           
                  "read_coverage_uniformity_score",
                  "junction_saturation_slope",
                  "median_insert_size",
                  "genes_contributing_to_80._of_reads",
                  "reads_mapping_sense_percentage",
                  "exonic_reads_minus_spike_ins",
                  "number_of_multimapped_reads",
                  "total_reads",
                  "number_of_uniquely_mapped_reads",
                  "mt_rna_pct", "mt_rrna_pct", "mt_trna_pct", "misc_rna_pct", "protein_coding_pct", "lncrna_pct", "snrna_pct", "snorna_pct", "spike_in_pct", "other_rna_biotypes_pct")

data <- read.delim("tables/sampleinfo_all-batches.tsv", header = TRUE, sep = "\t", fileEncoding = "UTF-8")

# Load metadata
metadata <- read.delim("tables/cfRNA-meta_per_sample_metadata.tsv", header = TRUE, sep = "\t", fill = TRUE)

# Strip decruyenaere samples of the "_*" in their sample_id
#data$sample_id[data$sequencing_batch == "decruyenaere"] <- 
  #sub("_.*", "", data$sample_id[data$sequencing_batch == "decruyenaere"])

# MOdify flomics_2 sample ids from FL-SAMP-ID TO SAMPID
#data$sequencing_batch[grepl("^FL", data$sequencing_batch)] <- "flomics_2"
#data$sample_id[data$sequencing_batch == "flomics_2"] <- 
  #gsub("_.*", "", data$sample_id[data$sequencing_batch == "flomics_2"])





# Step 1: Filter merged_df to keep only samples that appear in metadata
metadata_subset <- metadata[, c("run", "dataset_batch")]

filtered_df <- data[data$sample_id %in% metadata$run, ]
filtered_df <- merge(filtered_df, metadata_subset, by.x = "sample_id", by.y = "run", all.x = TRUE)


removed_samples <- data[!(data$sample_id %in% metadata$run), ]
print(removed_samples$sample_id)


cat("Original merged_df rows:", nrow(data), "\n") #should be 2458 
cat("Filtered to samples in metadata:", nrow(filtered_df), "\n") # should be 2360, if it's not, make sure you have NOT removed all Flomics_1 samples due to a mismatch between the metadata names and the sampleinfo from snakeda names :)



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
selected_columns <- c( "sample_id", "sample_name", "sequencing_batch", "status", "spike_in_pct", "protein_coding_pct", "percentage_of_spliced_reads", "dataset_batch.y", column_names)
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

table_filtered$log_genes_80 <- log(table_filtered$genes_contributing_to_80._of_reads)


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

bracket_df_2 <- data.frame(
  xmin = c("giraldez_phospho-rna-seq", "chalasani", "reggiardo_bioivt" , "block_150bp", "wei" ),
  xmax = c("wang", "toden", "reggiardo_dls" ,"rozowsky","wei" ),
  label = c("Custom", "Exome-based", "Whole RNA-Seq (oligo-dT pr.)", "Whole RNA-Seq (random pr.)", "cfDNA")
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
    
    offset <- 0.005
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


# ggplot_objects <- lapply(column_names, function(col_name) {
#   yvals <- table_filtered[[col_name]]
#   yvals <- yvals[!is.na(yvals) & is.finite(yvals)]
#   if (length(yvals) == 0) return(NULL)
#   
#   y_max <- max(yvals)
#   y_range <- diff(range(yvals))
#   bracket_y <- y_max + 0.05 * y_range
#   
#   bracket_df_top <- bracket_df %>%
#     mutate(
#       y.position = bracket_y,
#       y.position = ifelse(label == "Roskams-Hieter", bracket_y + 0.03 * y_range, y.position)
#     )
#   
#   label_block <- if (col_name == "mapped_percentage") {
#     labs(title = NULL, x = "Dataset", y = "% reads mapped to reference human genome")
#   } else if (col_name == "exonic_reads_minus_spike_ins") {
#     labs(title = NULL, x = "Dataset", y =  "Fraction of exonic reads (FER)")
#   } else if (col_name == "avg_mapped_read_length") {
#     labs(title = NULL, x = "Dataset", y =  "Effective fragment length\n(average mapped length)")
#   } else if (col_name == "spike_in_pct") {
#     labs(title = NULL, x = "Dataset", y =  "Fraction of reads mapping to ERCC spike-ins")
#   } else if (col_name == "percentage_of_spliced_reads") {
#     labs(title = NULL, x = "Dataset", y =  "Fraction of spliced reads (FSR)")
#   } else {
#     labs(title = NULL, x = "Dataset", y = clean_label(col_name))
#   }
#   
#   
#   p <- ggplot(table_filtered, aes(x = dataset_batch.y, y = .data[[col_name]], fill = dataset_batch.y)) +
#     geom_boxplot(aes(color = dataset_batch.y), alpha = 0.3, position = position_dodge(width = 0.75), outlier.shape = NA) +
#     geom_point(aes(y = .data[[col_name]], color = dataset_batch.y), 
#                position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.8), 
#                shape = 21, size = 1.5, stroke = 0.2, alpha = 0.6) +
#     label_block +
#     theme_classic() +
#     coord_cartesian(clip = "off") +
#     theme(text=element_text(family="Arial"),
#           panel.grid.major.x = element_blank(),
#               panel.grid.minor.x = element_blank(),
#               panel.grid.major.y = element_line(size = 0.8), 
#               panel.grid.minor.y = element_blank(),
#           plot.margin = margin(20, 20, 20, 50),
#       axis.text.x = element_text(angle = 45, hjust = 1, size = 10, vjust=0.85),
#           axis.title = element_text(size = 12, face = "bold"),
#           plot.title = element_blank(),
#           legend.position = "none") +
#     scale_x_discrete(labels = datasetsLabels) +
#     scale_fill_manual(values = datasetsPalette, labels = datasetsLabels) +
#     scale_color_manual(values = datasetsOutlinePalette, labels = datasetsLabels)
#   
#   p <- add_bottom_brackets(p, bracket_df, levels(table_filtered$dataset_batch.y), y_base = -0.055)
#   
#   p <- add_bottom_brackets(p, bracket_df_2, levels(table_filtered$dataset_batch.y), y_base = -0.025, col = "grey60", lwd=2)
#   
#   
#   return(p)
# })
# 
# 
# setwd("figures/")
# 
# for (i in 1:length(ggplot_objects)) {
#   col_name <- column_names[i]
#   output_file <- paste0(gsub(" ", "_", tolower(col_name)), "_external_datasets_boxplot_with_points.png")
#   ggsave(output_file, ggplot_objects[[i]], width = 11, height = 6, dpi = 600, device = ragg::agg_png)
#   output_file <- paste0(gsub(" ", "_", tolower(col_name)), "_external_datasets_boxplot_with_points.svg")
#   ggsave(output_file, ggplot_objects[[i]], width = 11, height = 6, dpi = 600, device = "svg")
# }

#######
# BPC prototype
#######

ggplot_objects <- lapply(column_names, function(col_name) {
  yvals <- table_filtered[[col_name]]
  yvals <- yvals[!is.na(yvals) & is.finite(yvals)]
  if (length(yvals) == 0) return(NULL)
  
  y_max <- max(yvals)
  y_range <- diff(range(yvals))
  bracket_y <- y_max + 0.05 * y_range
  
  bracket_df_top <- bracket_df %>%
    mutate(
      y.position = bracket_y,
      y.position = ifelse(label == "Roskams-Hieter", bracket_y + 0.03 * y_range, y.position)
    )
  
  label_block <- if (col_name == "mapped_percentage") {
    labs(title = NULL, x = "Dataset", y = "% reads mapped to reference human genome")
  } else if (col_name == "exonic_reads_minus_spike_ins") {
    labs(title = NULL, x = "Dataset", y =  "Fraction of exonic reads (FER)")
  } else if (col_name == "avg_mapped_read_length") {
    labs(title = NULL, x = "Dataset", y =  "Effective fragment length\n(average mapped length)")
  } else if (col_name == "spike_in_pct") {
    labs(title = NULL, x = "Dataset", y =  "Fraction of reads mapping to ERCC spike-ins")
  } else if (col_name == "percentage_of_spliced_reads") {
    labs(title = NULL, x = "Dataset", y =  "Fraction of spliced reads (FSR)")
  } else {
    labs(title = NULL, x = "Dataset", y = clean_label(col_name))
  }
  
  
  bpc_info <- metadata[, c("run", "broad_protocol_category")]
  table_filtered <- table_filtered %>%
    left_join(bpc_info, by = c("sample_id" = "run"))
  
  
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
  
  y_min <- min(table_filtered[[col_name]], na.rm = TRUE)
  y_max <- max(table_filtered[[col_name]], na.rm = TRUE)
  y_range <- y_max - y_min
  
  # Always place annotation tile below
  annotation_y <- 0 - 0.05 * y_range
  annotation_height <- 0.03 * y_range
  
  # Extend y-limits for tile 
  y_limit_min <- min(annotation_y - annotation_height, y_min - 0.02 * y_range)
  y_limit_max <- y_max + 0.1 * y_range
  
  
  
  bpc_order <- c("Custom", "Exome-based (EB)", "Whole RNA-Seq (oligo-dT pr.) (WRO)", "Whole RNA-Seq (random pr.) (WRR)", "cfDNA")
  
  annotation_df <- table_filtered %>%
    select(dataset_batch.y, broad_protocol_category) %>%
    distinct() %>%
    mutate(
      y = annotation_y,
      broad_protocol_category = factor(broad_protocol_category, levels = bpc_order)
    )
  
  
  
  p <- ggplot(table_filtered, aes(x = dataset_batch.y, y = .data[[col_name]], fill = dataset_batch.y)) +
    geom_boxplot(aes(color = dataset_batch.y), alpha = 0.3, position = position_dodge(width = 0.75), outlier.shape = NA) +
    geom_point(aes(y = .data[[col_name]], color = dataset_batch.y), 
               position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.8), 
               shape = 21, size = 1.5, stroke = 0.2, alpha = 0.6) +
    label_block +
    theme_classic() +
    coord_cartesian(clip = "off") +
    theme(text=element_text(family="Arial"),
          line = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(size = 0.8), 
          panel.grid.minor.y = element_blank(),
          plot.margin = margin(20, 20, 20, 50),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10, vjust=1.1),
          axis.title = element_text(size = 12, face = "bold"),
          axis.title.x = element_blank(),
          plot.title = element_blank(),
          legend.position = "none") +
    scale_x_discrete(labels = datasetsLabels) +
    scale_fill_manual(values = datasetsPalette, labels = datasetsLabels, guide = "none") +
    scale_color_manual(values = datasetsOutlinePalette, labels = datasetsLabels, guide= "none") +
    scale_y_continuous(
      limits = c(y_limit_min, y_limit_max)  # extend a bit for annotations
    )
  
  p <- p + new_scale_fill()  # VERY IMPORTANT, native ggplot does not like having two cals of scale_fill or scale_color, so ggnewscale is needed
  
  p <- p +
    geom_tile(data = annotation_df,
              aes(x = dataset_batch.y, y = y, fill = broad_protocol_category),
              width = 0.95, height = annotation_height,
              inherit.aes = FALSE) +
    scale_fill_manual(
      name = "Bropad Protocol Category (BPC)",
      values = bpc_colors,
      labels = bpc_labels
    ) + theme(text=element_text(family="Arial"),
              legend.position = "bottom",  
              legend.title = element_text(face = "bold", size = 9),
              legend.text = element_text(size = 9),
              legend.key.size = unit(0.4, "cm"),
              legend.spacing.x = unit(0.2, "cm"),
              legend.margin = margin(0, 0, 0, 0)
    )
  
  
  p <- add_bottom_brackets(p, bracket_df, levels(table_filtered$dataset_batch.y),  y_base = 0.03)
  #p <- add_bottom_brackets(p, bracket_df_2, levels(table_filtered$dataset_batch.y), y_base = -0.025, col = "grey60", lwd=2)
  
  return(p)
})

setwd("figures/")

for (i in 1:length(ggplot_objects)) {
  col_name <- column_names[i]
  output_file <- paste0(gsub(" ", "_", tolower(col_name)), "_annotation_row.png")
  ggsave(output_file, ggplot_objects[[i]], width = 11, height = 6, dpi = 600, device = ragg::agg_png)
  output_file <- paste0(gsub(" ", "_", tolower(col_name)), "_annotation_row.svg")
  ggsave(output_file, ggplot_objects[[i]], width = 11, height = 6, dpi = 600, device = "svg")
}
#####################################################
# Stats for manuscript SPLICED READS
####################################################
x <- table_filtered$percentage_of_spliced_reads[table_filtered$dataset_batch.y == "wei"]
median_value <- median(x, na.rm = TRUE)
sd_value <- sd(x, na.rm = TRUE)

cat("Median ± SD:", sprintf("%.2f ± %.2f", median_value, sd_value), "\n")


x <- table_filtered$percentage_of_spliced_reads[table_filtered$dataset_batch.y == "flomics_1"]
median_value <- median(x, na.rm = TRUE)
sd_value <- sd(x, na.rm = TRUE)

cat("Median ± SD:", sprintf("%.2f ± %.2f", median_value, sd_value), "\n")


x <- table_filtered$percentage_of_spliced_reads[table_filtered$dataset_batch.y == "flomics_2"]
median_value <- median(x, na.rm = TRUE)
sd_value <- sd(x, na.rm = TRUE)

cat("Median ± SD:", sprintf("%.2f ± %.2f", median_value, sd_value), "\n")

x <- table_filtered$percentage_of_spliced_reads[table_filtered$dataset_batch.y == "flomics_1"]
median_value <- median(x, na.rm = TRUE)
sd_value <- sd(x, na.rm = TRUE)

cat("Median ± SD:", sprintf("%.2f ± %.2f", median_value, sd_value), "\n")


x <- table_filtered$percentage_of_spliced_reads[table_filtered$dataset_batch.y == "flomics_2"]
median_value <- median(x, na.rm = TRUE)
sd_value <- sd(x, na.rm = TRUE)

cat("Median ± SD:", sprintf("%.2f ± %.2f", median_value, sd_value), "\n")

x <- table_filtered$percentage_of_spliced_reads[table_filtered$dataset_batch.y == "reggiardo_bioivt"]
median_value <- median(x, na.rm = TRUE)
sd_value <- sd(x, na.rm = TRUE)

cat("Median ± SD:", sprintf("%.2f ± %.2f", median_value, sd_value), "\n")


x <- table_filtered$percentage_of_spliced_reads[table_filtered$dataset_batch.y == "reggiardo_dls"]
median_value <- median(x, na.rm = TRUE)
sd_value <- sd(x, na.rm = TRUE)

cat("Median ± SD:", sprintf("%.2f ± %.2f", median_value, sd_value), "\n")

x <- table_filtered$percentage_of_spliced_reads[table_filtered$dataset_batch.y == "roskams_pilot"]
median_value <- median(x, na.rm = TRUE)
sd_value <- sd(x, na.rm = TRUE)

cat("Median ± SD:", sprintf("%.2f ± %.2f", median_value, sd_value), "\n")

x <- table_filtered$percentage_of_spliced_reads[table_filtered$dataset_batch.y == "roskams_validation"]
median_value <- median(x, na.rm = TRUE)
sd_value <- sd(x, na.rm = TRUE)

cat("Median ± SD:", sprintf("%.2f ± %.2f", median_value, sd_value), "\n")


#####################################################
# Stats for manuscript ng80
####################################################
x <- table_filtered$genes_contributing_to_80._of_reads[table_filtered$dataset_batch.y == "rozowsky"]
median_value <- median(x, na.rm = TRUE)
sd_value <- sd(x, na.rm = TRUE)

cat("Median ± SD:", sprintf("%.2f ± %.2f", median_value, sd_value), "\n")

x <- table_filtered$genes_contributing_to_80._of_reads[table_filtered$dataset_batch.y == "wei"]
median_value <- median(x, na.rm = TRUE)
sd_value <- sd(x, na.rm = TRUE)

cat("Median ± SD:", sprintf("%.2f ± %.2f", median_value, sd_value), "\n")

capture <- c("chalasani", "toden", "ibarra_buffy_coat", "ibarra_plasma_cancer", "ibarra_plasma_non_cancer", "ibarra_serum")

x <- table_filtered$genes_contributing_to_80._of_reads[table_filtered$dataset_batch.y %in% capture]
median_value <- median(x, na.rm = TRUE)
sd_value <- sd(x, na.rm = TRUE)

cat("Median ± SD:", sprintf("%.2f ± %.2f", median_value, sd_value), "\n")

gdna <- c("block_150bp", "block_300bp", "sun_2")

x <- table_filtered$genes_contributing_to_80._of_reads[table_filtered$dataset_batch.y == "sun_2"]
median_value <- median(x, na.rm = TRUE)
sd_value <- sd(x, na.rm = TRUE)

cat("Median ± SD:", sprintf("%.2f ± %.2f", median_value, sd_value), "\n")

gdna <- c("block_150bp", "block_300bp")

x <- table_filtered$genes_contributing_to_80._of_reads[table_filtered$dataset_batch.y %in% gdna]
median_value <- median(x, na.rm = TRUE)
sd_value <- sd(x, na.rm = TRUE)

cat("Median ± SD:", sprintf("%.2f ± %.2f", median_value, sd_value), "\n")

high_q <- c("chen", "decruyenaere", "flomics_1", "flomics_2", "moufarrej_site_1", "moufarrej_site_2", "reggiardo_bioivt", "reggiardo_dls", "roskams_pilot", "roskams_validation", "tao", "zhu")

x <- table_filtered$genes_contributing_to_80._of_reads[table_filtered$dataset_batch.y %in% high_q]
median_value <- median(x, na.rm = TRUE)
sd_value <- sd(x, na.rm = TRUE)

cat("Median ± SD:", sprintf("%.2f ± %.2f", median_value, sd_value), "\n")

#####################################################
# Stats for manuscript spike-ins
####################################################
summary(table_filtered$spike_in_pct[table_filtered$dataset_batch.y=="chen"])

######################################################
# High quality samples barplot
######################################################
# Create summary for high-quality sample barplot
setwd("~/fl-cfRNAmeta/")
quality_summary <- table_filtered %>%
  group_by(dataset_batch.y) %>%
  summarise(
    total = n(),
    high_quality = sum(genes_contributing_to_80._of_reads > 1000 & percentage_of_spliced_reads > 10, na.rm = TRUE)
  ) %>%
  mutate(
    percent_high_quality = 100 * high_quality / total,
    dataset_batch.y = factor(dataset_batch.y, levels = core_order)
  )

bpc_info <- metadata[, c("dataset_batch", "broad_protocol_category")]
bpc_info$dataset_batch.y <- bpc_info$dataset_batch


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

y_limit_min <- min(table_filtered[[col_name]], na.rm = TRUE) - 6
y_limit_max <- max(table_filtered[[col_name]], na.rm = TRUE) + 10

bpc_order <- c("Custom", "Exome-based (EB)", "Whole RNA-Seq (oligo-dT pr.) (WRO)", "Whole RNA-Seq (random pr.) (WRR)", "cfDNA")

annotation_df <- quality_summary %>%
  select(dataset_batch.y) %>%
  distinct() %>%
  left_join(bpc_info, by = "dataset_batch.y") %>%
  mutate(
    y = -3,  # row height below 0 (stacked bars start at 0)
    broad_protocol_category = factor(broad_protocol_category, levels = bpc_order)
  )


# Barplot of % high-quality samples per dataset batch
quality_plot <- ggplot(quality_summary, aes(x = dataset_batch.y, y = percent_high_quality, fill = dataset_batch.y)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_text(aes(label = paste0("N=", high_quality)), vjust = -0.5, size = 3) +
  labs(
    x = "Dataset",
    y = "Fraction of samples with NG80>1,000 and FSR>10%",
    title = "High-Quality Samples per Dataset"
  ) +
  #ylim(0, 105) +
  theme_classic() +
  coord_cartesian(clip = "off", ylim = c(-6, 105)) +
  theme(text=element_text(family="Arial"),
        line = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.8), 
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10, vjust=1.1),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank(),
        plot.title = element_blank(),
        plot.margin = margin(20, 20, 20, 45),
        legend = "none") +
  scale_x_discrete(labels = datasetsLabels) +
  scale_fill_manual(values = datasetsPalette, labels = datasetsLabels, guide = "none") +
  scale_color_manual(values = datasetsOutlinePalette, labels = datasetsLabels)

quality_plot <- quality_plot + new_scale_fill()

quality_plot <- quality_plot +
  geom_tile(
    data = annotation_df,
    aes(x = dataset_batch.y, y = y, fill = broad_protocol_category),
    width = 0.95, height = 3,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    name = "Broad Protocol Category (BPC)",
    values = bpc_colors,
    labels = bpc_labels
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.4, "cm"),
    legend.spacing.x = unit(0.2, "cm"),
    legend.margin = margin(0, 0, 0, 0)
  )


quality_plot <- add_bottom_brackets(quality_plot, bracket_df, levels(table_filtered$dataset_batch.y),  y_base = 0.03)
#quality_plot <- add_bottom_brackets(quality_plot, bracket_df_2, levels(table_filtered$dataset_batch.y), y_base = -0.025, col = "grey60", lwd=2)


# Save the new barplot
ggsave("figures/high_quality_sample_fraction_barplot.png", quality_plot, width = 11, height = 6, dpi = 600, device = ragg::agg_png)
ggsave("figures/high_quality_sample_fraction_barplot.svg", quality_plot, width = 11, height = 6, dpi = 600, device = "svg")


#######################################################
############### Fragments mapped to the expected strand
#######################################################

table_filtered$Fragments_mapping_to_expected_strand_pct <- 100 - as.numeric(table_filtered$reads_mapping_sense_percentage)

strandedness_info <- metadata[, c("run", "cdna_library_type")]
table_filtered <- table_filtered %>%
  left_join(strandedness_info, by = c("sample_id" = "run"))

# label_positions <- table_filtered %>%
#   group_by(dataset_batch.y, cdna_library_type) %>%
#   summarise(
#     y = max(Fragments_mapping_to_expected_strand_pct, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   mutate(y = y + 2)  


strandedness_colors <- c(
  "Reverse" = "#8DD3C7",
  "Unstranded" = "#BEBADA"
)

annotation_df <- table_filtered %>%
  select(dataset_batch.y, cdna_library_type) %>%
  distinct() %>%
  mutate(y = 105) 


# label_positions <- label_positions %>%
#   mutate(cdna_library_type = recode(cdna_library_type,
#                                     "Reverse" = "R",
#                                     "Unstranded" = "U"))



p <- ggplot(table_filtered, aes(x = dataset_batch.y, y = Fragments_mapping_to_expected_strand_pct, fill = dataset_batch.y)) +
  geom_boxplot(alpha = 0.3, aes(color = dataset_batch.y), position = position_dodge(width = 0.75), outlier.shape = NA) +
  geom_point(aes(y = Fragments_mapping_to_expected_strand_pct, color = dataset_batch.y), 
             position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.8), 
             shape = 21, size = 1.5, stroke = 0.2, alpha = 0.6) +
  #geom_hline(yintercept=100, linetype='dashed', col = 'lightgrey') +
  labs(title = "",
       x = "Dataset", y = "% fragments mapping to correct gene orientation") +
  theme_classic() +
  theme(text=element_text(family="Arial"),
        plot.margin = margin(20, 20, 20, 45),
    panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.8), 
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10, vjust = 0.85),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "none") +
  scale_x_discrete(labels = datasetsLabels) +
  scale_fill_manual(values = datasetsPalette, labels = datasetsLabels, guide = "none") +
  scale_color_manual(values = datasetsOutlinePalette, labels = datasetsLabels, guide = "none") +
  scale_y_continuous(
    limits = c(0, 110),  # extend a bit for annotations
    breaks = c(0, 25, 50, 75, 100),
    labels = function(x) paste0(x, "%")
  ) +
  coord_cartesian(clip = "off")

p <- p + new_scale_fill()  # VERY IMPORTANT, native ggplot does not like having two cals of scale_fill or scale_color, so ggnewscale is needed

p <- p +
  geom_tile(data = annotation_df,
            aes(x = dataset_batch.y, y = y, fill = cdna_library_type),
            width = 0.8, height = 4,
            inherit.aes = FALSE) +
  scale_fill_manual(
    name = "Library strandedness",
    values = strandedness_colors
  ) + theme(text=element_text(family="Arial"),
    legend.position = "top",  
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.4, "cm"),
    legend.spacing.x = unit(0.2, "cm"),
    legend.margin = margin(0, 0, 0, 0)
  )


p <- add_bottom_brackets(p, bracket_df, levels(table_filtered$dataset_batch.y),  y_base = -0.055)
p <- add_bottom_brackets(p, bracket_df_2, levels(table_filtered$dataset_batch.y), y_base = -0.025, col = "grey60", lwd=2)

ggsave("figures/fragments_mapped_expected_strand_with_strandedness_info.png", p, width = 11, height = 6, dpi = 600, device = ragg::agg_png)
ggsave("figures/fragments_mapped_expected_strand_with_strandedness_info.svg", p, width = 11, height = 6, dpi = 600, device = "svg")



###################################
################### Fragment number
###################################

table_filtered <- table_filtered %>%
  mutate(fragment_number = if_else(
    `dataset_batch.y` %in% c("giraldez_phospho-rna-seq", "giraldez_standard"),
    read_number,           # Single-end: keep as is
    read_number / 2        # Paired-end: divide by 2
  ))

p <- ggplot(table_filtered, aes(x = dataset_batch.y, y = fragment_number, fill = dataset_batch.y)) +
  geom_boxplot(alpha = 0.3, aes(color = dataset_batch.y), position = position_dodge(width = 0.75), outlier.shape = NA) +
  geom_point(aes(y = fragment_number, color = dataset_batch.y), 
             position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.8), 
             shape = 21, size = 1.5, stroke = 0.2, alpha = 0.6) +
  labs(title = "",
       x = "Dataset", y = "Fragment number") +
  theme_classic() +
  theme(text=element_text(family="Arial"),
        plot.margin = margin(20, 20, 20, 45),        
    panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank(),
                panel.grid.major.y = element_line(size = 0.8), 
                panel.grid.minor.y = element_blank(),
                axis.text.x = element_text(angle = 45, hjust = 1, size = 10, vjust = 0.85),
                axis.title = element_text(size = 12, face = "bold"),
                plot.title = element_text(size = 14, face = "bold"),
                legend.position = "none") +
  scale_x_discrete(labels = datasetsLabels) +
  scale_fill_manual(values = datasetsPalette, labels = datasetsLabels) +
  scale_color_manual(values = datasetsOutlinePalette, labels = datasetsLabels) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) + 
  coord_cartesian(clip = "off")

p <- add_bottom_brackets(p, bracket_df, levels(table_filtered$dataset_batch.y),  y_base = -0.055)
p <- add_bottom_brackets(p, bracket_df_2, levels(table_filtered$dataset_batch.y), y_base = -0.025, col = "grey60", lwd=2)

ggsave("figures/fragment_number.png", p, width = 11, height = 6, dpi = 600, device = ragg::agg_png)
ggsave("figures/fragment_number.svg", p, width = 11, height = 6, dpi = 600, device = "svg")


#######################################################
# NG80
#######################################################

p <- ggplot(table_filtered, aes(x = dataset_batch.y, y = log_genes_80, fill = dataset_batch.y)) +
  geom_boxplot(alpha = 0.3, aes(color = dataset_batch.y), position = position_dodge(width = 0.75), outlier.shape = NA) +
  geom_point(aes(y = log_genes_80, color = dataset_batch.y), 
             position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.8), 
             shape = 21, size = 1.5, stroke = 0.2, alpha = 0.6) +
  labs(title = "",
       x = "Dataset", y = "NG80") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "none") +
  scale_x_discrete(labels = datasetsLabels) +
  scale_fill_manual(values = datasetsPalette, labels = datasetsLabels) +
  scale_color_manual(values = datasetsOutlinePalette, labels = datasetsLabels)

y_breaks <- pretty(range(table_filtered$genes_contributing_to_80._of_reads, na.rm = TRUE), n = 10)
p <- p + scale_y_continuous(
  breaks = log(y_breaks),
  labels = round(y_breaks)
)

y_breaks <- log(c(100, 500, 1000, 5000, 10000, 20000))

y_labels <- c(100, 500, 1000, 5000, 10000, 20000)

p <- ggplot(table_filtered, aes(x = dataset_batch.y, y = log_genes_80, fill = dataset_batch.y)) +
  geom_boxplot(alpha = 0.3, aes(color = dataset_batch.y), position = position_dodge(width = 0.75), outlier.shape = NA) +
  geom_point(aes(y = log_genes_80, color = dataset_batch.y), 
             position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.8), 
             shape = 21, size = 1.5, stroke = 0.2, alpha = 0.6) +
  labs(title = "",
       x = "Dataset", y = "NG80") +
  theme_classic() +
  theme(text=element_text(family="Arial"),
        plot.margin = margin(20, 20, 20, 45),      
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.8), 
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10, vjust = 0.85),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "none") +
  scale_x_discrete(labels = datasetsLabels) +
  scale_fill_manual(values = datasetsPalette, labels = datasetsLabels) +
  scale_color_manual(values = datasetsOutlinePalette, labels = datasetsLabels) +
  scale_y_continuous(breaks = y_breaks, labels = y_labels) + 
  coord_cartesian(clip = "off")

p <- add_bottom_brackets(p, bracket_df, levels(table_filtered$dataset_batch.y),  y_base = -0.055)
p <- add_bottom_brackets(p, bracket_df_2, levels(table_filtered$dataset_batch.y), y_base = -0.025, col = "grey60", lwd=2)


ggsave("figures/ng80_non_transformed_axis.png", p, width = 11, height = 6, dpi = 600, device = ragg::agg_png)
ggsave("figures/ng80_non_transformed_axis.svg", p, width = 11, height = 6, dpi = 600, device = "svg")

#################################
# NG80 protein coding
################################

ng_only_mrna <- read.delim("tables/genes_contributing_to_percentage_reads.tsv")

ng80_table <- table_filtered %>%
  left_join(ng_only_mrna, by = c("sample_name" = "Sample"))

ng80_table$log_genes_80_pc <- log(ng80_table$number_of_genes_contributing_to_80._of_reads)

y_breaks <- log(c(100, 500, 1000, 5000, 10000, 20000))

y_labels <- c(100, 500, 1000, 5000, 10000, 20000)

p <- ggplot(ng80_table, aes(x = dataset_batch.y, y = log_genes_80, fill = dataset_batch.y)) +
  geom_boxplot(alpha = 0.3, aes(color = dataset_batch.y), position = position_dodge(width = 0.75), outlier.shape = NA) +
  geom_point(aes(y = log_genes_80, color = dataset_batch.y), 
             position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.8), 
             shape = 21, size = 1.5, stroke = 0.2, alpha = 0.6) +
  labs(title = "",
       x = "Dataset", y = "NG80 (protein coding genes)") +
  theme_classic() +
  theme(text=element_text(family="Arial"),
        plot.margin = margin(20, 20, 20, 45),  
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.8), 
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10, vjust = 0.85),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "none") +
  scale_x_discrete(labels = datasetsLabels) +
  scale_fill_manual(values = datasetsPalette, labels = datasetsLabels) +
  scale_color_manual(values = datasetsOutlinePalette, labels = datasetsLabels) +
  scale_y_continuous(breaks = y_breaks, labels = y_labels) + 
  coord_cartesian(clip = "off")

p <- add_bottom_brackets(p, bracket_df, levels(table_filtered$dataset_batch.y),  y_base = -0.055)
p <- add_bottom_brackets(p, bracket_df_2, levels(table_filtered$dataset_batch.y), y_base = -0.025, col = "grey60", lwd=2)


ggsave("figures/ng80_mrna_non_transformed_axis.png", p, width = 11, height = 6, dpi = 600, device = ragg::agg_png)
ggsave("figures/ng80_mrna_non_transformed_axis.svg", p, width = 11, height = 6, dpi = 600, device = "svg")


################################
# NpcG80 vs NG80
###############################

ng80_table$ratio <- ng80_table$number_of_genes_contributing_to_80._of_reads / ng80_table$genes_contributing_to_80._of_reads

y_breaks <- log(c(100, 500, 1000, 5000, 10000, 20000))

y_labels <- c(100, 500, 1000, 5000, 10000, 20000)

to_keep <- c("chalasani", "toden", "ibarra_buffy_coat", "ibarra_plasma_cancer", "ibarra_plasma_non_cancer", "ibarra_serum", "chen", "decruyenaere", "flomics_1", "flomics_2", "moufarrej_site_1", "moufarrej_site_2", "reggiardo_bioivt", "reggiardo_dls", "roskams_pilot", "roskams_validation", "tao", "zhu")

core_order_filtered <- c("chalasani", "ibarra_buffy_coat", "ibarra_plasma_cancer", "ibarra_plasma_non_cancer", "ibarra_serum", "toden", "reggiardo_bioivt", "reggiardo_dls", "chen", "decruyenaere", "flomics_1", "flomics_2", "moufarrej_site_1", "moufarrej_site_2", "roskams_pilot", "roskams_validation", "tao", "zhu")               


ng80_table <- ng80_table %>%
  filter(dataset_batch.y %in% to_keep) %>%
  mutate(dataset_batch.y = fct_relevel(factor(dataset_batch.y), core_order_filtered ))

bracket_df <- data.frame(
  xmin = c("ibarra_buffy_coat", "moufarrej_site_1", "reggiardo_bioivt", "roskams_pilot"),
  xmax = c( "ibarra_serum","moufarrej_site_2", "reggiardo_dls", "roskams_validation"),
  label = c( "Ibarra", "Moufarrej", "Reggiardo", "Roskams-Hieter")
)

bracket_df_2 <- data.frame(
  xmin = c("chalasani", "reggiardo_bioivt" , "chen"),
  xmax = c("toden", "reggiardo_dls" ,"zhu"),
  label = c("Exome-based", "Whole RNA-Seq (oligo-dT pr.)", "Whole RNA-Seq (random pr.)")
)

p <- ggplot(ng80_table, aes(x = dataset_batch.y, y = ratio, fill = dataset_batch.y)) +
  geom_boxplot(alpha = 0.3, aes(color = dataset_batch.y), position = position_dodge(width = 0.75), outlier.shape = NA) +
  geom_point(aes(y = ratio, color = dataset_batch.y), 
             position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.8), 
             shape = 21, size = 1.5, stroke = 0.2, alpha = 0.6) +
  labs(title = "",
       x = "Dataset", y = "NP80 / NG80") +
  theme_classic() +
  theme(text=element_text(family="Arial"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.8), 
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10, vjust = 0.85),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "none") +
  ylim(c(0,2.5))+
  geom_hline(yintercept=1, linetype='dashed', col = 'darkgrey') +
  scale_x_discrete(labels = datasetsLabels) +
  scale_fill_manual(values = datasetsPalette, labels = datasetsLabels) +
  scale_color_manual(values = datasetsOutlinePalette, labels = datasetsLabels) +
  #scale_y_continuous(breaks = y_breaks, labels = y_labels) + 
  coord_cartesian(clip = "off")

add_bottom_brackets_filtered <- function(p, bracket_df, x_levels, offset = 0.005, y_base = -0.03, height = 0.015, col = "black", lwd = 0.95) {
  for (i in seq_len(nrow(bracket_df))) {
    x1 <- which(x_levels == bracket_df$xmin[i])
    x2 <- which(x_levels == bracket_df$xmax[i])
    if (length(x1) == 0 || length(x2) == 0) next
    
    # Add small spacing between brackets
    offset <- 0.005 #spacing
    x_start <- ((x1 - 1) / length(core_order_filtered)) + offset  # LEFT tick
    x_end   <- (x2 / length(core_order_filtered)) - offset  
    
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


p <- add_bottom_brackets_filtered(p, bracket_df, levels(ng80_table$dataset_batch.y),  y_base = -0.055)
p <- add_bottom_brackets_filtered(p, bracket_df_2, levels(ng80_table$dataset_batch.y), y_base = -0.025, col = "grey60", lwd=2)


ggsave("figures/ng80_ratio_non_transformed_axis_filtered.png", p, width = 11, height = 6, dpi = 600, device = ragg::agg_png)
ggsave("figures/ng80_ratio_non_transformed_axis_filtered.svg", p, width = 11, height = 6, dpi = 600, device = "svg")


capture <- c("chalasani", "toden", "ibarra_buffy_coat", "ibarra_plasma_cancer", "ibarra_plasma_non_cancer", "ibarra_serum")

x <- ng80_table$ratio[ng80_table$dataset_batch.y %in% capture]
median_value <- median(x, na.rm = TRUE)
sd_value <- sd(x, na.rm = TRUE)

cat("Median ± SD:", sprintf("%.2f ± %.2f", median_value, sd_value), "\n")

high_q <- c("chen", "decruyenaere", "flomics_1", "flomics_2", "moufarrej_site_1", "moufarrej_site_2", "reggiardo_bioivt", "reggiardo_dls", "roskams_pilot", "roskams_validation", "tao", "zhu")

x <- ng80_table$ratio[ng80_table$dataset_batch.y %in% high_q]
median_value <- median(x, na.rm = TRUE)
sd_value <- sd(x, na.rm = TRUE)

cat("Median ± SD:", sprintf("%.2f ± %.2f", median_value, sd_value), "\n")
#################################
########### Diversity scatterplot
#################################
p <- ggplot(data = table_filtered, aes(x = percentage_of_spliced_reads, y = genes_contributing_to_80._of_reads, color = dataset_batch.y)) +
  geom_point(size = 2, alpha = 0.7) + 
  geom_density_2d(aes(color = dataset_batch.y), alpha = 0.5, size = 0.8) +
  labs(title = "",
       x = "Percentage of spliced reads",
       y = "NG80",
       color = "")+
  theme_minimal(base_size = 14) + 
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5), 
    axis.title = element_text(size = 14, face = "bold"), 
    axis.text = element_text(size = 12), 
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10), 
    legend.position = "right", 
    legend.key.height = unit(0.5, "lines"),  
    legend.box = "vertical",
    panel.grid.major = element_line(color = "gray80"), 
    panel.grid.minor = element_blank(),
    plot.background = element_rect(
      fill = "white",
      colour = "white"
    )
  ) +
  scale_color_manual(values = datasetsPalette, labels= datasetsLabels) 

ps <- p  +  scale_y_continuous(trans=log10_trans()) +
  guides(color = guide_legend(ncol = 1))


ggsave("figures/diversity_scatterplot.png", plot = ps, width = 12, height = 8, dpi = 600, device = ragg::agg_png)
ggsave("figures/diversity_scatterplot.svg", plot = ps, width = 12, height = 8,  dpi = 600, device = "svg")


##################################
# Per dataset diversity scatterplot
##################################
adjusted_palette <- datasetsPalette

table_filtered$dataset_batch.y <- factor(table_filtered$dataset_batch.y, levels = core_order)
adjusted_palette <- adjusted_palette[core_order]



# Set proper order and labels
table_filtered$dataset_batch.y <- factor(
  table_filtered$dataset_batch.y,
  levels = core_order,
  labels = datasetsLabels[core_order]
)

adjusted_palette <- adjusted_palette[core_order]
names(adjusted_palette) <- datasetsLabels[core_order]  # to match new factor labels




p_facet <- ggplot(
  data = table_filtered,
  aes(x = percentage_of_spliced_reads,
      y = genes_contributing_to_80._of_reads,
      color = dataset_batch.y)) +
  geom_point(size = 2, alpha = 0.8) +
  facet_wrap(~ dataset_batch.y, scales = "fixed") +
  scale_color_manual(
    values = adjusted_palette,
    drop = FALSE
  ) +
  labs(
    x = "Fraction of spliced reads (FSR)",
    y = "NG80",
    color = "Dataset"
  ) +
  theme_minimal(base_size = 13) +
  theme(text=element_text(family="Arial"),
    strip.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.text = element_text(size = 9),
    legend.title = element_text(face = "bold"),
    legend.key.height = unit(0.5, "lines"),
    plot.background = element_rect(fill = "white", colour = "white")
  ) +
  guides(color = guide_legend(ncol = 1)) +
  scale_y_continuous(trans = log10_trans())


# Define axis limits globally
x_lim <- range(table_filtered$percentage_of_spliced_reads, na.rm = TRUE)
y_lim <- range(table_filtered$genes_contributing_to_80._of_reads, na.rm = TRUE)

p_facet <- ggplot(
  data = table_filtered,
  aes(x = percentage_of_spliced_reads,
      y = genes_contributing_to_80._of_reads,
      color = dataset_batch.y)) +
  geom_point(size = 2, alpha = 0.8) +
  # Trick to enforce shared limits even with scales = "free"
  geom_blank(aes(x = x_lim[1], y = y_lim[1])) +
  geom_blank(aes(x = x_lim[2], y = y_lim[2])) +
  facet_wrap(~ dataset_batch.y, scales = "free") +
  scale_color_manual(values = adjusted_palette, drop = FALSE) +
  labs(
    x = "Fraction of spliced reads (FSR)",
    y = "NG80",
    color = "Dataset"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.text = element_text(size = 9),
    legend.title = element_text(face = "bold"),
    legend.key.height = unit(0.5, "lines"),
    plot.background = element_rect(fill = "white", colour = "white"),
    panel.spacing = unit(0.8, "lines")
  ) +
  guides(color = guide_legend(ncol = 1)) +
  scale_y_continuous(trans = log10_trans())


ggsave("figures/diversity_scatterplot_facet.png", p_facet, width = 20, height = 10, dpi = 600, device = ragg::agg_png)
ggsave("figures/diversity_scatterplot_facet.svg", p_facet, width = 20, height = 10, device = "svg")


#################################
# Per dataset k2 results vs mapping rate
################################
k2_results <- read.delim("~/fl-cfRNAmeta/tables/taxa_simple_df_w_batch.tsv")

# First ensure dataset_batch is available in k2_results
# Then clean sample_name only for affected dataset_batch values
# k2_results_clean <- k2_results %>%
#   mutate(sample_name = case_when(
#     dataset_batch %in% c("flomics_2", "decruyenaere") ~ sub("_.*", "", sample_name),
#     TRUE ~ sample_name
#   ))

# Compute microbial reads again if not already there
k2_results$microbial <- k2_results$fungi + k2_results$bacteria

# Merge into table_filtered using cleaned sample names
table_filtered <- table_filtered %>%
  left_join(k2_results %>% select(sample_name, microbial),
            by = c("sample_id" = "sample_name"))


p_microbial <- ggplot(
  data = table_filtered,
  aes(x = mapped_percentage,
      y = microbial,
      color = dataset_batch.y)) +
  
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, color = "grey", alpha = 0.8) +
  
  stat_cor(
    method = "pearson",
    label.x.npc = "left", 
    label.y.npc = "top",   
    size = 3,
    aes(label = ..r.label..),
    color = "black"
  ) +
  
  facet_wrap(~ dataset_batch.y, scales = "free") +
  
  scale_color_manual(
    values = adjusted_palette,
    drop = FALSE
  ) +
  labs(
    x = "% reads mapped to reference human genome",
    y = "% microbial reads",
    color = "Dataset"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    text=element_text(family="Arial"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(size = 0.8), 
    panel.grid.minor.y = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.text = element_text(size = 9),
    legend.title = element_text(face = "bold"),
    legend.key.height = unit(0.5, "lines"),
    plot.background = element_rect(fill = "white", colour = "white")
  ) +
  guides(color = guide_legend(ncol = 1))


ggsave("figures/microbial_vs_mapped.png", p_microbial, width = 20, height = 10, dpi = 600, device = ragg::agg_png)
ggsave("figures/microbial_vs_mapped.svg", p_microbial, width = 20, height = 10, device = "svg")


####################################
############ Biotype stacked barplot (DEPRECATED)
####################################
biotype_columns <- grep("_pct$", names(table_filtered), value = TRUE)

biotype_data <- merged_df[, 66:147]
biotype_data <- biotype_data %>% dplyr::select(-contains('_fc'))

# Create supergroups with corrected names
biotype_data$protein_coding_1 <- rowSums(biotype_data[, c("ig_c_gene", "ig_d_gene", "ig_j_gene",
                                                          "ig_v_gene", "tr_c_gene", "tr_d_gene", "tr_j_gene",
                                                          "tr_v_gene", "protein_coding")], na.rm = TRUE)

#biotype_data$mt_ncrna <- rowSums(biotype_data[, c("mt_rrna", "mt_trna")], na.rm = TRUE)

biotype_data$smrna <- rowSums(biotype_data[, c("mirna", "snrna", "snorna")], na.rm = TRUE)

biotype_data$pseudogene_1 <- rowSums(biotype_data[, c(
  "ig_c_pseudogene", "ig_j_pseudogene", "ig_v_pseudogene", "ig_pseudogene", "tr_j_pseudogene",
  "tr_v_pseudogene", "polymorphic_pseudogene", "processed_pseudogene", "pseudogene", "rrna_pseudogene",
  "transcribed_processed_pseudogene", "transcribed_unitary_pseudogene", "transcribed_unprocessed_pseudogene",
  "translated_processed_pseudogene", "translated_unprocessed_pseudogene", "unitary_pseudogene",
  "unprocessed_pseudogene")], na.rm = TRUE)


biotype_data$other <- rowSums(biotype_data[, c("tec", "ribozyme", "srna", "scrna", "scarna", "vault_rna")], na.rm = TRUE)

# Select supergroup columns + sequencing_batch
biotype_data_selected <- biotype_data %>%
  mutate(sequencing_batch = merged_df$sequencing_batch) %>%
  dplyr::select(sequencing_batch,
                protein_coding_1,
                lncrna,
                mt_rrna,
                mt_trna,
                smrna,
                rrna,
                misc_rna,
                pseudogene_1,
                spike_in,
                other)

biotype_data_selected$sequencing_batch[grepl("^FL", biotype_data_selected$sequencing_batch)] <- "liquidx"


# Summarize across datasets
biotype_summary_supergroups <- biotype_data_selected %>%
  group_by(sequencing_batch) %>%
  summarise(across(everything(), ~mean(.x, na.rm = TRUE))) %>%
  pivot_longer(-sequencing_batch, names_to = "biotype", values_to = "mean_pct")


# Clean labels
biotype_summary_supergroups$biotype <- gsub("_1$", "", biotype_summary_supergroups$biotype)
biotype_summary_supergroups$biotype <- gsub("_", " ", biotype_summary_supergroups$biotype)

# Use your dataset label mapping
biotype_summary_supergroups$sequencing_batch <- factor(biotype_summary_supergroups$sequencing_batch,
                                                       levels = names(datasetsLabels),
                                                       labels = datasetsLabels)

biotype_summary_supergroups$sequencing_batch <- factor(biotype_summary_supergroups$sequencing_batch,
                                                       levels = core_order)


# Distinct palette for supergroups
supergroup_palette <- c(
  "protein coding" = "#1f77b4",
  "lncrna"         = "#ff7f0e",
  "mt rrna"        = "#2ca02c",
  "mt trna"        = "darkgreen",
  "smrna"          = "#d62728",
  "rrna"           = "#9467bd",
  "misc rna"       = "#17becf",
  "pseudogene"     = "#bcbd22",
  "spike in"       = "#8c564b",
  "other"          = "gray"
)

# Plot
p_supergroups <- ggplot(biotype_summary_supergroups, aes(x = sequencing_batch, y = mean_pct, fill = biotype)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Dataset", y = "% of reads", fill = "Biotype Group") +
  scale_fill_manual(values = supergroup_palette) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        legend.key.size = unit(0.5, "cm"),
        legend.spacing.x = unit(0.3, 'cm'))

ggsave("biotype_supergroup_distribution_per_dataset.png", p_supergroups, width = 12, height = 4, dpi = 600)
ggsave("biotype_supergroup_distribution_per_dataset.pdf", p_supergroups, width = 12, height = 4, dpi = 600)
########################################

###############################
############ Protein coding pct DEPRECATED
###############################
p <- ggplot(table_filtered, aes(x = sequencing_batch, y = protein_coding_pct, fill = sequencing_batch)) +
  geom_boxplot(alpha = 0.6, aes(color = sequencing_batch), position = position_dodge(width = 0.75), outlier.shape = NA) +
  geom_point(aes(y = protein_coding_pct, color = sequencing_batch), 
             position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.8), 
             shape = 21, size = 1.5, stroke = 0.2, alpha = 0.6) +
  labs(title = "",
       x = "Dataset", y = "Protein coding %") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "none") +
  scale_x_discrete(labels = datasetsLabels) +
  scale_fill_manual(values = datasetsPalette, labels = datasetsLabels) +
  scale_color_manual(values = datasetsOutlinePalette, labels = datasetsLabels)


ggsave("protein_coding_pct.png", p, width = 9, height = 6, dpi = 600)
ggsave("protein_coding_pct.pdf", p, width = 9, height = 6, dpi = 600)

# p <- ggplot(data = table_filtered, aes(x = Exonic_percentage, y = log(Genes_contributing_to_80._of_reads), color = dataset)) +
#   geom_point() +
#   geom_density_2d() + # comment to remove the density lines
#   labs(x = "exonic percentage",
#        y = "log(# of genes contributing to 80% of reads)",
#        color = "") +
#   theme_minimal() +
#   scale_color_manual(values = datasetsPalette, labels = datasetsLabels)+
#   scale_fill_manual(values = datasetsPalette) 

##################################
################## Shannon entropy DEPRECATED
##################################

shannon <- read.table("~/cfRNA-meta/exp_mat/shannon_entropy_with_flomics_1_andliquidx.csv", header = TRUE, sep = ",")
colnames(shannon) <- c("sequencing_batch", "sample_id", "ShannonEntropy")

table_filtered_clean <- table_filtered
table_filtered_clean$sample_id <- gsub("FL-|-", "", table_filtered_clean$sample_id)


table_filtered_combined <- left_join(table_filtered_clean, shannon, by = "sample_id")

p <- ggplot(data = table_filtered_combined, aes(x = exonic_reads_minus_spike_ins, y = ShannonEntropy, color = dataset)) +
  geom_point(size = 2, alpha = 0.7) + 
  geom_density_2d(aes(color = dataset), alpha = 0.5, size = 0.8) +
  labs(title = "",
       x = "Exonic percentage minus spike-ins",
       y = "Shannon Entropy",
       color = "")+
  theme_minimal(base_size = 14) + 
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5), 
    axis.title = element_text(size = 14, face = "bold"), 
    axis.text = element_text(size = 12), 
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10), 
    legend.position = "right", 
    panel.grid.major = element_line(color = "gray80"), 
    panel.grid.minor = element_blank(),
    plot.background = element_rect(
      fill = "white",
      colour = "white"
    )
  ) +
  scale_color_manual(values = datasetsPalette, labels= datasetsLabels) 

p <- ggplot(table_filtered_combined, aes(x = sequencing_batch.x, y = ShannonEntropy, fill = sequencing_batch.x)) +
  geom_boxplot(alpha = 0.3, aes(color = sequencing_batch.x), position = position_dodge(width = 0.75), outlier.shape = NA ) +
  geom_point(aes(y = ShannonEntropy, color = sequencing_batch.x), 
             position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.8), 
             shape = 16, size = 2) +
  labs(title = "",
       x = "Dataset", y = "Shannon Entropy") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        legend.spacing.y = unit(0.2, 'cm')) +
  scale_x_discrete(labels = datasetsLabels) +
  scale_fill_manual(values = datasetsPalette, labels = datasetsLabels) +
  scale_color_manual(values = datasetsOutlinePalette, labels = datasetsLabels)

ggsave("shannon_entropy_boxplot.png", plot = p, width = 9, height = 6, dpi = 600)
ggsave("shannon_entropy_boxplot.pdf", plot = p, width = 9, height = 6, dpi = 600)



################################################################################
# FOR SLIDE DECK
################################################################################
table_filtered_slide_deck <- table_filtered

# Remove datasets
table_filtered_slide_deck <- table_filtered_slide_deck[table_filtered_slide_deck$sequencing_batch != "chalasani", ]
table_filtered_slide_deck <- table_filtered_slide_deck[table_filtered_slide_deck$sequencing_batch != "toden", ]
table_filtered_slide_deck <- table_filtered_slide_deck[table_filtered_slide_deck$sequencing_batch != "ibarra", ]
table_filtered_slide_deck <- table_filtered_slide_deck[table_filtered_slide_deck$sequencing_batch != "sun_1", ]
table_filtered_slide_deck <- table_filtered_slide_deck[table_filtered_slide_deck$sequencing_batch != "sun_2", ]
table_filtered_slide_deck <- table_filtered_slide_deck[table_filtered_slide_deck$sequencing_batch != "reggiardo", ]
table_filtered_slide_deck <- table_filtered_slide_deck[table_filtered_slide_deck$sequencing_batch != "wang", ]
table_filtered_slide_deck <- table_filtered_slide_deck[table_filtered_slide_deck$sequencing_batch != "giraldez", ]
table_filtered_slide_deck <- table_filtered_slide_deck[table_filtered_slide_deck$sequencing_batch != "moufarrej", ]
table_filtered_slide_deck <- table_filtered_slide_deck[table_filtered_slide_deck$sequencing_batch != "rozowsky", ]

#joint_table_final <- joint_table_final[joint_table_final$sequencing_batch != "taowei", ]

table_filtered_slide_deck$sequencing_batch <- as.character(table_filtered_slide_deck$sequencing_batch)

table_filtered_slide_deck$dataset <- ifelse(
  table_filtered_slide_deck$sequencing_batch %in% c("block_1", "block_2"), "block",
  ifelse(table_filtered_slide_deck$sequencing_batch %in% c("roskams_1", "roskams_2"), "roskams",
         table_filtered_slide_deck$sequencing_batch)
)

table_filtered_slide_deck$dataset <- factor(table_filtered_slide_deck$dataset, levels = c( "liquidx","block" ,"decruyenaere", "zhu", "chen",  "ngo", "roskams", "moufarrej","tao","rozowsky", "taowei"))

datasetsPalette=c( "liquidx" = "#144d6b", "block" = "#b3b3b3" , "decruyenaere" =  "#b3b3b3" , "zhu" = "#b3b3b3", "chen" = "#b3b3b3", "ngo" =  "#b3b3b3", "roskams" = "#b3b3b3" , "moufarrej" = "#b3b3b3",  "tao" ="#b3b3b3","rozowsky" = "#006600", "taowei"="#B32400")

datasetsLabels=c("liquidx" = "Flomics", "block" = "Block", "chen" = "Chen", "decruyenaere" = "Decruyenaere", "ngo" = "Ngo", "roskams" = "Roskams-Hieter","moufarrej" = "Moufarrej ","rozowsky"="ENCODE\n(bulk tissue RNA-Seq)", "tao" = "Tao", "zhu" ="Zhu", "taowei"="Wei (cfDNA)")


darken_color <- function(color, factor = 1.3) {
  rgb_col <- col2rgb(color) / 255
  darker_rgb <- pmin(rgb_col * (1 / factor), 1)
  rgb(darker_rgb[1], darker_rgb[2], darker_rgb[3])
}

datasetsOutlinePalette <- sapply(datasetsPalette, darken_color)

p <- ggplot(table_filtered_slide_deck, aes(x = dataset, y = percentage_of_spliced_reads, fill = dataset)) +
  geom_boxplot(aes(color = dataset), alpha = 0.6, position = position_dodge(width = 0.75), outlier.shape = NA) +
  geom_point(aes(y = percentage_of_spliced_reads, color = dataset), 
             position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.8), 
             shape = 21, size = 1.5, stroke = 0.2, alpha = 0.6) +
  labs(title = "",
       x = "Dataset", y = "Percent of spliced reads") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "none", 
        legend.title = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        legend.spacing.y = unit(0.2, 'cm')) +
  scale_x_discrete(labels = datasetsLabels) +
  scale_fill_manual(values = datasetsPalette, labels = datasetsLabels) +
  scale_color_manual(values = datasetsOutlinePalette, labels = datasetsLabels) + # Match outline colors to fill
  guides(fill = guide_legend(override.aes = list(color = "black"), ncol = 1)) 

ggsave("percent_of_spliced_reads_slide_deck_darker_outline.png", p, width = 8, height = 6)
ggsave("percent_of_spliced_reads_slide_deck_darker_outline.pdf", p, width = 8, height = 6)

darken_color <- function(color, factor = 1.3) {
  rgb_col <- col2rgb(color) / 255
  darker_rgb <- pmin(rgb_col * (1 / factor), 1)
  rgb(darker_rgb[1], darker_rgb[2], darker_rgb[3])
}

datasetsOutlinePalette <- sapply(datasetsPalette, darken_color)


p <- ggplot(table_filtered_slide_deck, aes(x = dataset, y = genes_contributing_to_80._of_reads, fill = dataset)) +
  geom_boxplot(aes(color = dataset), alpha = 0.6, position = position_dodge(width = 0.75), outlier.shape = NA) +
  geom_point(aes(y = genes_contributing_to_80._of_reads, color = dataset), 
             position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.8), 
             shape = 21, size = 1.5, stroke = 0.2, alpha = 0.6) +
  labs(title = "",
       x = "Dataset", y = "Library diversity\n(# genes contributing to 80% of reads)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "none", 
        legend.title = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        legend.spacing.y = unit(0.2, 'cm')) +
  scale_x_discrete(labels = datasetsLabels) +
  scale_fill_manual(values = datasetsPalette, labels = datasetsLabels) +
  scale_color_manual(values = datasetsOutlinePalette, labels = datasetsLabels) + # Match outline colors to fill
  guides(fill = guide_legend(override.aes = list(color = "black"), ncol = 1)) 

ggsave("library_diversity_slide_deck.png", p, width = 8, height = 6)
ggsave("library_diversity_slide_deck.pdf", p, width = 8, height = 6)

ps <- p  +  scale_y_continuous(trans=log10_trans()) 
ggsave("library_diversity_slide_deck_log10.png", ps, width = 8, height = 6)
ggsave("library_diversity_slide_deck_log10.pdf", ps, width = 8, height = 6)


################## Shannon entropy

shannon <- read.table("~/cfRNA-meta/exp_mat/shannon_entropy_with_liquidx.csv", header = TRUE, sep = ",")
colnames(shannon) <- c("sequencing_batch", "sample_id", "ShannonEntropy")

# Remove FL- and - from sample_id column so both tables match
table_filtered_slide_deck$sample_id <- gsub("-", "", gsub("^FL-", "", table_filtered_slide_deck$sample_id))

table_filtered_combined <- left_join(table_filtered_slide_deck, shannon, by = "sample_id")

p <- ggplot(data = table_filtered_combined, aes(x = percentage_of_spliced_reads, y = genes_contributing_to_80._of_reads, color = dataset)) +
  geom_point(size = 2, alpha = 0.7) + 
  geom_density_2d(aes(color = dataset), alpha = 0.5, size = 0.8) +
  labs(title = "",
       x = "Percent of spliced reads",
       y = "Library diversity",
       color = "")+
  theme_minimal(base_size = 14) + 
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5), 
    axis.title = element_text(size = 14, face = "bold"), 
    axis.text = element_text(size = 12), 
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10), 
    legend.position = "right", 
    panel.grid.major = element_line(color = "gray80"), 
    panel.grid.minor = element_blank(),
    plot.background = element_rect(
      fill = "white",
      colour = "white"
    )
  ) +
  scale_color_manual(values = datasetsPalette, labels= datasetsLabels) 
ps <- p  +  scale_y_continuous(trans=log2_trans()) 

ggsave("library_vs_spliced.png", plot = ps, width = 12, height = 8, dpi = 600)
ggsave("library_vs_spliced.pdf", plot = ps, width = 12, height = 8, dpi = 600)


p <- ggplot(data = table_filtered_combined, aes(x = percentage_of_spliced_reads, y = genes_contributing_to_80._of_reads, color = dataset)) +
  geom_point(size = 2, alpha = 0.7) + 
  geom_density_2d(aes(color = dataset), alpha = 0.5, size = 0.8) +
  labs(title = "",
       x = "Percent of spliced reads",
       y = "Library diversity",
       color = "") +
  theme_minimal(base_size = 14) + 
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5), 
    axis.title = element_text(size = 14, face = "bold"), 
    axis.text = element_text(size = 12), 
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10), 
    legend.position = "left", 
    panel.grid.major = element_line(color = "gray80"), 
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", colour = "white")
  ) +
  scale_color_manual(values = datasetsPalette, labels = datasetsLabels) +
  scale_y_continuous(trans = log2_trans())

pd <- ggMarginal(p, groupColour = TRUE, groupFill = TRUE, type = "density")
pb <- ggMarginal(p, groupColour = FALSE, groupFill = TRUE, type = "boxplot")

ggsave("library_vs_spliced_bp.png", plot = pb, width = 15, height = 8, dpi = 600)
ggsave("library_vs_spliced_bp.pdf", plot = pb, width = 15, height = 8, dpi = 600)


datasetsLabels <- c(
  "liquidx" = "Flomics",
  "block" = "Block",
  "decruyenaere" = "Decruyenaere",
  "zhu" = "Zhu",
  "chen" = "Chen",
  "ngo" = "Ngo",
  "roskams" = "Roskams-Hieter",
  "moufarrej" = "Moufarrej",
  "tao" = "Tao",
  "rozowsky" = "ENCODE (bulk tissue RNA-Seq)",
  "taowei" = "Wei (cfDNA)"
)
datasetsPalette=c( "Flomics" = "#144d6b",
                   "Block" = "#b3b3b3" ,
                   "Decruyenaere" =  "#b3b3b3",
                   "Zhu" = "#b3b3b3",
                   "Chen" = "#b3b3b3",
                   "Ngo" =  "#b3b3b3", 
                   "Roskams-Hieter" = "#b3b3b3",
                   "Moufarrej" = "#b3b3b3",  
                   "Tao" ="#b3b3b3",
                   "ENCODE (bulk tissue RNA-Seq)" = "#006600", 
                   "Wei (cfDNA)"="#B32400")



table_plot <- table_filtered_combined %>%
  mutate(
    dataset = datasetsLabels[dataset],  # relabel dataset values
    log_lib_diversity = genes_contributing_to_80._of_reads
    #log_lib_diversity = log10(genes_contributing_to_80._of_reads)
  ) %>%
  filter(!is.na(percentage_of_spliced_reads), !is.na(log_lib_diversity))

order <- c("Flomics", "Block", "Decruyenaere", "Zhu", "Chen", "Ngo", "Roskams-Hieter", "Moufarrej", "Tao", "Wei (cfDNA)", "ENCODE (bulk tissue RNA-Seq)")
table_plot$dataset <- factor(table_plot$dataset, levels = order)


p <- ggplot(table_plot, aes(x = percentage_of_spliced_reads, y = log_lib_diversity, color = dataset)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_density_2d(alpha = 0.5, size = 0.8) +
  
  geom_xsideboxplot(
    aes(y = dataset, fill = dataset), orientation = "y",
    width = 0.6,
    outlier.shape = NA,
    #coef = 0,
    alpha = 0.3,
    side = "top",
    show.legend = FALSE
  ) +
  
  geom_ysideboxplot(
    aes(x = dataset, group = dataset, fill = dataset), orientation = "x",
    width = 0.6,
    outlier.shape = NA,
    #coef = 0,
    alpha = 0.3,
    side = "right",
    show.legend = FALSE
  ) +
  
  scale_color_manual(values = datasetsPalette, labels = datasetsLabels) +
  scale_fill_manual(values = datasetsPalette) +
  
  labs(
    x = "Percent of spliced reads",
    y = "Library diversity\n(# of genes contributing to 80% of reads)",
    color = ""
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    ggside.panel.scale = 0.3, 
    ggside.axis.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    ggside.axis.text.x = element_blank(),
    ggside.panel.scale.x = 0.25,  # smaller panel
    ggside.panel.scale.y = 0.25, 
    plot.background = element_rect(fill = "white", colour = "white"),
    panel.grid.major = element_blank()
  )

ps <- p + scale_y_continuous(trans=log10_trans()) 

ggsave("library_vs_spliced_boxplots.png", plot = ps, width = 13, height = 10, dpi = 600)
ggsave("library_vs_spliced_boxplots.pdf", plot = ps, width = 13, height = 10, dpi = 600)


p <- ggplot(table_filtered_combined, aes(x = dataset, y = ShannonEntropy, fill = dataset)) +
  geom_boxplot(alpha = 0.3, color = "black", position = position_dodge(width = 0.75), outlier.shape = NA ) +
  geom_point(aes(y = ShannonEntropy, color = dataset), 
             position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.8), 
             shape = 16, size = 2) +
  labs(title = "",
       x = "Dataset", y = "Shannon Entropy") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        legend.spacing.y = unit(0.2, 'cm')) +
  scale_x_discrete(labels = datasetsLabels) +
  scale_fill_manual(values = datasetsPalette, labels = datasetsLabels) +
  scale_color_manual(values = datasetsPalette, labels = datasetsLabels)

ggsave("shannon_entropy_boxplot.png", plot = p, width = 9, height = 6, dpi = 600)
ggsave("shannon_entropy_boxplot.pdf", plot = p, width = 9, height = 6, dpi = 600)

