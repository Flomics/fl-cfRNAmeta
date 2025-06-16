library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(RColorBrewer)
library(ComplexHeatmap)
library(purrr)
library(jsonlite)
library(showtext)
font_add("DejaVu Sans", regular = "DejaVuSans.ttf")
showtext_opts(dpi = 600)  # MUST come before showtext_auto()
showtext_auto()
theme_set(theme_classic(base_family = "DejaVu Sans"))

data_heatmap <- read.table("tables/cfRNA-meta_per_batch_metadata.tsv", header = TRUE, sep = "\t", na.strings = c("", "NA"))

data_heatmap$centrifugation_step_1 <- as.character(data_heatmap$centrifugation_step_1)
data_heatmap$centrifugation_step_2 <- as.character(data_heatmap$centrifugation_step_2)
data_heatmap$centrifugation_step_1 <- trimws(data_heatmap$centrifugation_step_1)
data_heatmap$centrifugation_step_2 <- trimws(data_heatmap$centrifugation_step_2)


# Drop n_samples column and the GC1 and GC5, and other unnecessary columns
data_heatmap <- subset(data_heatmap, select = -c(n_samples, libraryselection, genes_contributing_to_1._of_reads, genes_contributing_to_5._of_reads, cancer_stage,
                                                 phenotype_subtype, stage_cancer_simple, is_excluded_from_study, race, lab, dev_stage, description, assemblyname, resequenced_sample,
                                                 project, organism, bioproject, datastore_provider, datastore_region))

giraldez_batches <- c("giraldez_standard", "giraldez_phospho-rna-seq")

data_heatmap$library_prep_kit_short_name[data_heatmap$dataset_batch %in% giraldez_batches] <- 
  data_heatmap$library_prep_kit_short[data_heatmap$dataset_batch %in% giraldez_batches]


metadata_long <- data_heatmap %>%
  pivot_longer(-dataset_batch, names_to = "variable", values_to = "value") %>%
  mutate(value = ifelse(is.na(value), "NA", value))  # Replace NA with string

metadata_long$variable <- factor(metadata_long$variable)

metadata_matrix <- metadata_long %>%
  pivot_wider(names_from = dataset_batch, values_from = value) %>%
  column_to_rownames("variable")

n_values <- apply(metadata_matrix, 1, function(x) length(unique(x)))


get_palette_with_na <- function(varname, base) {
  values <- unique(as.character(metadata_matrix[varname, ]))
  if (varname %in% c("centrifugation_step_1", "centrifugation_step_2")) {
    desired_order <- c("1000g", "1500g", "1600g", "1900g", "1940g", "2000g", "2500g", "3000g", 
                       "3400g", "6000g", "12000g", "13000g", "15000g", "16000g",
                       "Unspecified", "placeholder", "None")
    values <- intersect(desired_order, values)
    real_values <- values[!values %in% c("NA", "Unspecified", "placeholder")]
  } else if (varname == "read_length") {
    desired_order <- c("1x50", "1x75", "2x75", "2x100", "2x150", "NA")
    values <- intersect(desired_order, unique(as.character(metadata_matrix[varname, ])))
    real_values <- values[values != "NA"]
  } else {
    values <- sort(unique(as.character(metadata_matrix[varname, ])))
    real_values <- values[values != "NA"]
  }
  
  
  real_values <- values[values != "NA"]
  n_real <- length(real_values)
  
  if (!(base %in% rownames(brewer.pal.info))) {
    stop(paste("Invalid RColorBrewer palette:", base))
  }
  
  max_colors <- brewer.pal.info[base, "maxcolors"]
  if (n_real == 0) {
    return(c("NA" = "grey80"))
  }
  
  n_use <- min(n_real, max_colors)
  color_fun <- colorRampPalette(brewer.pal(n_use, base))
  palette <- color_fun(n_real)
  
  color_map <- setNames(palette, real_values)
  color_map["NA"] <- "grey80"
  color_map["Unspecified"] <- "grey60"
  color_map["placeholder"] <- "lavenderblush1"
  color_map["None"] <- "grey70"
  
  return(color_map)
}


mappings <- fromJSON("src/dataset_mappings.json")

clean_dataset_names <- unlist(mappings$datasetsLabels)
core_order <- unlist(mappings$datasetVisualOrder)


# all names alphabetically except "rozowsky" and "wei" last
ordered_datasets <- c(clean_dataset_names[core_order])


# Rename columns to clean display names
colnames(metadata_matrix) <- clean_dataset_names[colnames(metadata_matrix)]


metadata_matrix <- metadata_matrix[, ordered_datasets]


metadata_matrix <- apply(metadata_matrix, c(1, 2), function(x) {
  x <- trimws(x)
  if (tolower(x) == "unspecified") "Unspecified" else x
}) %>% as.data.frame(check.names = FALSE)

palette_list <- list(
  biomaterial = get_palette_with_na("biomaterial", "Set1"),
  nucleic_acid_type = get_palette_with_na("nucleic_acid_type", "PiYG"),
  library_selection = get_palette_with_na("library_selection", "Paired"),
  rna_extraction_kit_short_name = get_palette_with_na("rna_extraction_kit_short_name", "Set3"),
  library_prep_kit_short_name = get_palette_with_na("library_prep_kit_short_name", "Paired"),
  dnase = get_palette_with_na("dnase", "Accent"),
  cdna_library_type  = get_palette_with_na("cdna_library_type", "Set3"),
  read_length = get_palette_with_na("read_length", "Greens"),
  "centrifugation_step_1" = get_palette_with_na("centrifugation_step_1", "Reds") ,
  "centrifugation_step_2" = get_palette_with_na("centrifugation_step_2", "Oranges")
)


clean_names <- c(
  biomaterial = "Biomaterial",
  nucleic_acid_type = "Nucleic acid type",
  library_selection = "Library selection",
  rna_extraction_kit_short_name = "RNA extraction kit",
  library_prep_kit_short_name = "Library prep kit",
  dnase = "DNAse treatment",
  cdna_library_type = "cDNA library type",
  read_length = "Read length",
  centrifugation_step_1 = "Centrifugation, step 1",
  centrifugation_step_2 = "Centrifugation, step 2"
)

row_order <- c(
  "centrifugation_step_1",
  "centrifugation_step_2",
  "biomaterial",
  "nucleic_acid_type",
  "rna_extraction_kit_short_name",
  "dnase",
  "library_prep_kit_short_name",
  "library_selection",
  "cdna_library_type",
  "read_length"
)

heatmap_list <- lapply(row_order, function(var) {
  if (is.null(palette_list[[var]])) {
    message("Skipping ", var, ": no palette defined.")
    return(NULL)
  }
  
  values <- as.character(metadata_matrix[var, ])
  
  # Replace actual NA with "NA" string
  values[is.na(values)] <- "NA"
  
  if (var %in% c("centrifugation_step_1", "centrifugation_step_2")) {
    values <- factor(values, levels  = c("1000g", "1500g", "1600g", "1900g", "1940g", "2000g", "2500g", "3000g", 
                      "3400g", "6000g", "12000g", "13000g", "15000g", "16000g",
                      "Unspecified", "placeholder", "None"))
    values <- as.character(values)
  }
  
  
  if (var == "read_length") {
    values <- factor(values, levels = c("1x50", "1x75", "2x75", "2x100", "2x150", "NA"))
    values <- as.character(values)
  }
  
  unique_vals <- unique(values)
  color_vector <- palette_list[[var]]
  
  if (is.null(color_vector) || is.function(color_vector)) {
    stop(paste("No valid palette for variable:", var))
  }
  color_map <- color_vector[unique_vals]
  
  names(color_map) <- unique_vals
  
  if (!"NA" %in% names(color_map)) {
    color_map["NA"] <- "grey80"
  }
  color_map["Unspecified"] <- "grey60"
  color_map["placeholder"] <- "lavenderblush1"
  color_map["None"] <- "grey70"
  if (var == "read_length") {
    legend_order <- c("1x50", "1x75", "2x75", "2x100", "2x150", "NA")
    color_map <- color_map[intersect(legend_order, names(color_map))]
  }
  if (var %in% c("centrifugation_step_1", "centrifugation_step_2")) {
    legend_order <- c("1000g", "1500g", "1600g", "1900g", "1940g", "2000g", "2500g", "3000g", 
                      "3400g", "6000g", "12000g", "13000g", "15000g", "16000g",
                      "Unspecified", "placeholder", "None")
    color_map <- color_map[intersect(legend_order, names(color_map))]
  }
  
  # Always move special values to the end of the legend
  special_levels <- c("Unspecified", "placeholder", "NA", "None")
  ordered_vals <- c(
    setdiff(names(color_map), special_levels),
    intersect(special_levels, names(color_map))
  )
  color_map <- color_map[ordered_vals]
  
  var_pretty <- clean_names[[var]]
  
  ordered_pretty_names <- clean_names[row_order]

  mat <- matrix(values, nrow = 1, dimnames = list(factor(var_pretty, levels = ordered_pretty_names), colnames(metadata_matrix)))

  
  Heatmap(
    mat,
    name = var_pretty,
    col = color_map,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_names_rot = 45,
    row_names_side = "left",
    height = unit(1, "cm"),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.rect(
        x = x, y = y,
        width = width, height = height,
        gp = gpar(fill = fill, col = "white", lwd = 0.5)
      )
    }
    
  )
}) %>% discard(is.null)

names(heatmap_list) <- row_order


ht_list <- Reduce(`%v%`, heatmap_list)


ragg::agg_png("figures/fig_1b_metadata_heatmap.png", width = 18, height = 10, units = "in", res = 600)
draw(ht_list, heatmap_legend_side = "right")
dev.off()

pdf("figures/fig_1b_metadata_heatmap.pdf", width = 18, height = 10)
showtext::showtext_begin()
draw(ht_list, heatmap_legend_side = "right")
showtext::showtext_end()
dev.off()



