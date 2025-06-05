library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(RColorBrewer)
library(ComplexHeatmap)
library(purrr)

data_heatmap <- read.table("../tables/cfRNA-meta_per_batch_metadata.tsv", header = TRUE, sep = "\t", na.strings = c("", "NA"))

# Drop n_samples column
data_heatmap <- subset(data_heatmap, select = -c(n_samples))

# Pivot + replace NA with "NA" string before anything else
metadata_long <- data_heatmap %>%
  pivot_longer(-dataset_batch, names_to = "variable", values_to = "value") %>%
  mutate(value = ifelse(is.na(value), "NA", value))  # Replace NA with string

metadata_long$variable <- factor(metadata_long$variable)

metadata_matrix <- metadata_long %>%
  pivot_wider(names_from = dataset_batch, values_from = value) %>%
  column_to_rownames("variable")

# Recalculate number of levels per variable (including "NA")
n_values <- apply(metadata_matrix, 1, function(x) length(unique(x)))


get_palette_with_na <- function(varname, base) {
  values <- unique(as.character(metadata_matrix[varname, ]))
  if (varname == "read_length") {
    desired_order <- c("1x50", "1x75", "2x75", "2x100", "2x150", "NA")
    values <- intersect(desired_order, c(values, "NA"))  # preserve actual values only, in desired order
  } else {
    values <- sort(values)
  }
  
  real_values <- values[values != "NA"]
  n_real <- length(real_values)
  
  # Defensive: Check palette availability
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
  return(color_map)
}


palette_list <- list(
  biomaterial = get_palette_with_na("biomaterial", "Set1"),
  nucleic_acid_type = get_palette_with_na("nucleic_acid_type", "PiYG"),
  libraryselection = get_palette_with_na("libraryselection", "Paired"),
  rna_extraction_kit_short_name = get_palette_with_na("rna_extraction_kit_short_name", "Set2"),
  library_prep_kit_short_name = get_palette_with_na("library_prep_kit_short_name", "Dark2"),
  dnase = get_palette_with_na("dnase", "Accent"),
  cdna_library_type  = get_palette_with_na("cdna_library_type", "Set3"),
  read_length = get_palette_with_na("read_length", "Greens")
)


clean_names <- c(
  biomaterial = "Biomaterial",
  nucleic_acid_type = "Nucleic acid type",
  libraryselection = "Library selection",
  rna_extraction_kit_short_name = "RNA extraction kit",
  library_prep_kit_short_name = "Library prep kit",
  dnase = "DNAse treatment",
  cdna_library_type = "cDNA library type",
  read_length = "Read length"
)

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
  reggiardo = "Reggiardo"
)

colnames(metadata_matrix) <- clean_dataset_names[colnames(metadata_matrix)]

heatmap_list <- lapply(rownames(metadata_matrix), function(var) {
  if (is.null(palette_list[[var]])) {
    message("Skipping ", var, ": no palette defined.")
    return(NULL)
  }
  
  values <- as.character(metadata_matrix[var, ])
  
  # Replace actual NA with "NA" string
  values[is.na(values)] <- "NA"
  
  # Build color palette including grey for "NA"
  unique_vals <- unique(values)
  palette <- palette_list[[var]]
  
  # Expand the palette with grey for "NA"
  color_map <- palette[unique_vals]
  names(color_map) <- unique_vals
  
  if (!"NA" %in% names(color_map)) {
    color_map["NA"] <- "grey80"
  }
  
  var_pretty <- clean_names[[var]]
  
  # Build matrix
  mat <- matrix(values, nrow = 1, dimnames = list(var_pretty, colnames(metadata_matrix)))
  
  
  # Create heatmap
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

draw(ht_list, heatmap_legend_side = "right")

