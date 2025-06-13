library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(RColorBrewer)
library(ComplexHeatmap)
library(purrr)

data_heatmap <- read.table("tables/cfRNA-meta_per_batch_metadata.tsv", header = TRUE, sep = "\t", na.strings = c("", "NA"))

# Drop n_samples column and the GC1 and GC5, and other unnecessary columns
data_heatmap <- subset(data_heatmap, select = -c(n_samples, libraryselection, genes_contributing_to_1._of_reads, genes_contributing_to_5._of_reads))

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
  if (varname == "read_length") {
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
  return(color_map)
}


clean_dataset_names <- c(
  chen = "Chen",
  zhu = "Zhu",
  roskams_pilot = "Roskams-Hieter (pilot)",
  roskams_validation = "Roskams-Hieter (validation)",
  ngo = "Ngo",
  ibarra_serum = "Ibarra (serum)",
  ibarra_plasma = "Ibarra (plasma)",
  ibarra_buffy_coat = "Ibarra (buffy coat)",
  toden = "Toden",
  chalasani = "Chalasani",
  block_150bp = "Block (150bp)",
  block_300bp = "Block (300bp)",
  rozowsky = "ENCODE\n(bulk tissue RNA-Seq)",
  tao = "Tao",
  wei = "Wei (cfDNA)",
  moufarrej = "Moufarrej",
  wang = "Wang",
  giraldez_standard = "Giráldez (standard)",
  "giraldez_phospho-rna-seq" = "Giráldez (phospho-RNA-seq)",
  sun_2 = "Sun",
  decruyenaere = "Decruyenaere",
  reggiardo = "Reggiardo",
  flomics_1 = "Flomics 1",
  flomics_2 = "Flomics 2"
)
# Rename columns to clean display names
colnames(metadata_matrix) <- clean_dataset_names[colnames(metadata_matrix)]

core_order <- setdiff(names(clean_dataset_names), c("rozowsky", "wei"))
ordered_datasets <- c(sort(clean_dataset_names[core_order]), clean_dataset_names[c("rozowsky", "wei")])

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
  read_length = get_palette_with_na("read_length", "Greens")
)


clean_names <- c(
  biomaterial = "Biomaterial",
  nucleic_acid_type = "Nucleic acid type",
  library_selection = "Library selection",
  rna_extraction_kit_short_name = "RNA extraction kit",
  library_prep_kit_short_name = "Library prep kit",
  dnase = "DNAse treatment",
  cdna_library_type = "cDNA library type",
  read_length = "Read length"
)


heatmap_list <- lapply(rownames(metadata_matrix), function(var) {
  if (is.null(palette_list[[var]])) {
    message("Skipping ", var, ": no palette defined.")
    return(NULL)
  }
  
  values <- as.character(metadata_matrix[var, ])
  
  # Replace actual NA with "NA" string
  values[is.na(values)] <- "NA"
  
  if (var == "read_length") {
    values <- factor(values, levels = c("1x50", "1x75", "2x75", "2x100", "2x150", "NA"))
    values <- as.character(values)
  }
  
  unique_vals <- unique(values)
  palette <- palette_list[[var]]
  
  color_map <- palette[unique_vals]
  names(color_map) <- unique_vals
  
  if (!"NA" %in% names(color_map)) {
    color_map["NA"] <- "grey80"
  }
  color_map["Unspecified"] <- "grey60"
  if (var == "read_length") {
    legend_order <- c("1x50", "1x75", "2x75", "2x100", "2x150", "NA")
    color_map <- color_map[intersect(legend_order, names(color_map))]
  }
  
  var_pretty <- clean_names[[var]]
  
  mat <- matrix(values, nrow = 1, dimnames = list(var_pretty, colnames(metadata_matrix)))
  
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


ht_list <- Reduce(`%v%`, heatmap_list)


png("figures/fig_1b_metadata_heatmap.png", width = 18, height = 10, res = 600, units = "in", bg = "white")
draw(ht_list, heatmap_legend_side = "right")
dev.off()

