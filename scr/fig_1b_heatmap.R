library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(RColorBrewer)
library(ComplexHeatmap)

data_heatmap <- read.table("../tables/fig_1b.tsv", header = TRUE, sep = "\t")

metadata_long <- data_heatmap %>%
  pivot_longer(-Study.name, names_to = "variable", values_to = "value")

metadata_long$variable <- factor(metadata_long$variable)

metadata_matrix <- metadata_long %>%
  pivot_wider(names_from = Study.name, values_from = value) %>%
  column_to_rownames("variable")

get_palette <- function(n, base = "Set3") {
  colorRampPalette(brewer.pal(min(12, brewer.pal.info[base, "maxcolors"]), base))(n)
}

n_values <- apply(metadata_matrix, 1, function(x) length(unique(na.omit(x))))

palette_list <- list(
  Biomaterial = get_palette(n_values["Biomaterial"], base = "Set1"),
  Nucleic.acid.type = get_palette(n_values["Nucleic.acid.type"], base = "Pastel1"),
  RNA.extraction.kit = get_palette(n_values["RNA.extraction.kit"], base = "Set2"),
  Library.prep.kit = get_palette(n_values["Library.prep.kit"], base = "Dark2"),
  DNAse.treatment = get_palette(n_values["DNAse.treatment"], base = "Accent"),
  Strandedness = get_palette(n_values["Strandedness"], base = "Set3"),
  Library.layout = get_palette(n_values["Library.layout"], base = "Paired"),
  Read.length = get_palette(n_values["Read.length"], base = "Spectral")
)


heatmap_list <- lapply(rownames(metadata_matrix), function(var) {
  values <- as.character(metadata_matrix[var, ])
  palette <- palette_list[[var]]
  
  unique_vals <- unique(na.omit(values))
  color_map <- setNames(palette[seq_along(unique_vals)], unique_vals)
  
  mat <- matrix(values, nrow = 1, dimnames = list(var, colnames(metadata_matrix)))
  
  Heatmap(
    mat,
    name = var,
    col = color_map,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_names_rot = 45,
    row_names_side = "left",
    height = unit(1, "cm")
  )
})

ht_list <- Reduce(`%v%`, heatmap_list)

draw(ht_list, heatmap_legend_side = "right")


