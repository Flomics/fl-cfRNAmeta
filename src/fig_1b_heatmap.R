library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(RColorBrewer)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
library(purrr)
library(jsonlite)
library(wesanderson)
library(grid)
library(showtext)
library(svglite)
library("extrafont")
loadfonts()

setwd("~/fl-cfRNAmeta/")

data_heatmap <- read.table("tables/cfRNA-meta_per_batch_metadata.tsv", header = TRUE, sep = "\t")
data_heatmap <- data_heatmap %>%
  mutate(
    dnase = case_when(
      dnase == "" ~ "None",
      TRUE ~ dnase
    ),
    rna_extraction_kit_short_name = case_when(
      rna_extraction_kit_short_name == "" ~ "None",
      TRUE ~ rna_extraction_kit_short_name
    )
  )



data_heatmap$centrifugation_step_1 <- as.character(data_heatmap$centrifugation_step_1)
data_heatmap$centrifugation_step_2 <- as.character(data_heatmap$centrifugation_step_2)
data_heatmap$centrifugation_step_1 <- trimws(data_heatmap$centrifugation_step_1)
data_heatmap$centrifugation_step_2 <- trimws(data_heatmap$centrifugation_step_2)


# Drop n_samples column and the GC1 and GC5, and other unnecessary columns
data_heatmap <- subset(data_heatmap, select = -c(n_samples, libraryselection, cancer_stage,
                                                 phenotype_subtype, stage_cancer_simple, is_excluded_from_study, race, lab, dev_stage, description, assemblyname,
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

label_replacements <- c(
  "PAXgene blood ccfDNA" = "PAXgene blood\nccfDNA",
  "Maxwell (plasma and serum)" = "Maxwell (plasma and\nserum)",
  "Qiagen miRNeasy Plasma" = "Qiagen miRNeasy\nPlasma",
  "Qiagen RNeasy Plasma" = "Qiagen RNeasy\nPlasma",
  "Qiagen RNeasy Tissue" = "Qiagen RNeasy\nTissue",
  "DNase I + HL-dsDNase" = "DNase I +\nHL-dsDNase",
  "Illumina TruSeq small RNA" = "Illumina TruSeq small\nRNA",
  "PNK-treated Illumina TruSeq small RNA" = "PNK-treated Illumina\nTruSeq small RNA",
  "Whole-Exome Capture" = "Whole-Exome\nCapture",
  "Whole-Transcriptome < 80 nt" = "Whole-Transcriptome\n< 80 nt",
  "Whole-Transcriptome > 200 nt" = "Whole-Transcriptome\n> 200 nt",
  "Whole RNA-Seq (oligo-dT pr.) (WRO)" = "Whole RNA-Seq\n(oligo-dT pr.) (WRO)",
  "Whole RNA-Seq (random pr.) (WRR)" = "Whole RNA-Seq\n(random pr.) (WRR)"
)

for (old_label in names(label_replacements)) {
  metadata_matrix[metadata_matrix == old_label] <- label_replacements[[old_label]]
}


n_values <- apply(metadata_matrix, 1, function(x) length(unique(x)))



get_palette_with_na <- function(varname, base, expand = TRUE) {
  values <- unique(as.character(metadata_matrix[varname, ]))
  
  if (varname %in% c("centrifugation_step_1", "centrifugation_step_2")) {
    desired_order <- c("1000g", "1500g", "1600g", "1900g", "1940g", "2000g", "2500g", "3000g", 
                       "3400g", "6000g", "12000g", "13000g", "15000g", "16000g",
                       "Unspecified", "placeholder", "None")
    values <- intersect(desired_order, values)
  } else if (varname == "read_length") {
    desired_order <- c("1x50", "1x75", "2x75", "2x100", "2x150", "NA")
    values <- intersect(desired_order, values)
  } else {
    values <- sort(values)
  }
  
  real_values <- values[!values %in% c("NA", "Unspecified", "placeholder", "None")]
  n_real <- length(real_values)
  
  if (n_real == 0) return(c("NA" = "grey80"))
  
  
  # Special handling for centrifugation palettes
  else if (varname %in% c("centrifugation_step_1", "centrifugation_step_2")) {
    palette <- colorRampPalette(brewer.pal(9, base)[3:9])(n_real)
  }
  
  # Special handling for read lenght
  else if (varname %in% c("read_length")) {
    palette <- colorRampPalette(brewer.pal(9, base)[3:9])(n_real)
  }
  # Default behavior
  else {
    if (base %in% names(wesanderson::wes_palettes)) {
      base_colors <- wesanderson::wes_palettes[[base]]
      if (n_real <= length(base_colors)) {
        palette <- wesanderson::wes_palette(base, n_real, type = "discrete")
      } else if (expand) {
        palette <- colorRampPalette(base_colors)(n_real)
      } else {
        stop(paste("wes_palette", base, "only supports", length(base_colors), "colors, but", n_real, "requested. Set expand=TRUE to interpolate."))
      }
    } else {
      if (!(base %in% rownames(RColorBrewer::brewer.pal.info))) {
        stop(paste("Invalid palette:", base))
      }
      max_colors <- RColorBrewer::brewer.pal.info[base, "maxcolors"]
      n_use <- min(n_real, max_colors)
      palette <- colorRampPalette(RColorBrewer::brewer.pal(n_use, base))(n_real)
    }
  }
  
  color_map <- setNames(palette, real_values)
  
  # Add fallback/special values
  color_map["NA"] <- "grey80"
  color_map["Unspecified"] <- "grey20"
  color_map["placeholder"] <- "lavenderblush1"
  color_map["None"] <- "grey70"
  
  return(color_map)
}

mappings <- fromJSON("src/dataset_mappings.json")

clean_dataset_names <- unlist(mappings$datasetsLabels)
core_order <- unlist(mappings$datasetVisualOrder)

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
  rna_extraction_kit_short_name = get_palette_with_na("rna_extraction_kit_short_name", "Spectral", expand = TRUE),
  library_prep_kit_short_name = get_palette_with_na("library_prep_kit_short_name", "Paired"),
  dnase = get_palette_with_na("dnase", "Accent"),
  cdna_library_type  = get_palette_with_na("cdna_library_type", "Set3"),
  read_length = get_palette_with_na("read_length", "Greens"),
  "centrifugation_step_1" = get_palette_with_na("centrifugation_step_1", "Reds") ,
  "centrifugation_step_2" = get_palette_with_na("centrifugation_step_2", "Oranges"),
  plasma_tubes_short_name = get_palette_with_na("plasma_tubes_short_name", "GrandBudapest2"),
  broad_protocol_category = get_palette_with_na("broad_protocol_category", "Cavalcanti1")
)

# I could not include into the function the special handling of smarter kits, changing their colours manually
palette_list$library_prep_kit_short_name[6] <- "#FDBE85"
palette_list$library_prep_kit_short_name[7] <- "#FD8D3C"
palette_list$library_prep_kit_short_name[8] <- "#E6550D"

palette_list$broad_protocol_category[1] <- "#AECAD9"
palette_list$broad_protocol_category[2] <- "#D9BBAE"
palette_list$broad_protocol_category[3] <- "#BBE0BB"
palette_list$broad_protocol_category[4] <- "#D0AED9"
palette_list$broad_protocol_category[5] <- "#D9D6AE"

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
  centrifugation_step_2 = "Centrifugation, step 2",
  plasma_tubes_short_name = "Blood collection tube",
  broad_protocol_category = "Broad protocol\ncategory (BPC)"
)

row_order <- c(
  "plasma_tubes_short_name",
  "centrifugation_step_1",
  "centrifugation_step_2",
  "biomaterial",
  "nucleic_acid_type",
  "rna_extraction_kit_short_name",
  "dnase",
  "library_prep_kit_short_name",
  "library_selection",
  "cdna_library_type",
  "read_length",
  "broad_protocol_category"
)

#Hacky way to transform Reggiardo dataset into a single batch, in terms of protocol
metadata_matrix$`Reggiardo (BioIVT)` <- NULL
colnames(metadata_matrix)[11] <- "Reggiardo"

bracket_df <- data.frame(
  xmin = c("Block (2x75bp)", "Giráldez (phospho)", "Ibarra (buffy coat)", "Moufarrej (Site 1)", "Roskams (pilot)"),
  xmax = c("Block (2x150bp)", "Giráldez (standard)", "Ibarra (serum)", "Moufarrej (Site 2)", "Roskams (validation)"),
  label = c("Block", "Giráldez", "Ibarra", "Moufarrej", "Roskams")
)

core_order <- colnames(metadata_matrix)  # your x-axis order

# Convert dataset names to indices
bracket_df$xmin_idx <- match(bracket_df$xmin, core_order) 
bracket_df$xmax_idx <- match(bracket_df$xmax, core_order) 

library(circlize)

# Define fallback colors globally so we can reuse them
fallback_colors <- c(
  "Unspecified" = "grey20",
  #"placeholder" = "lavenderblush1",
  "None" = "grey70"#,
  #"NA" = "grey80"
)

# Create custom legend for centrifugation
make_centrifugation_legend1 <- function(title, range = c(1000, 12000), colors = c("#FCBBA1", "#67000D")) {
  col_fun <- colorRamp2(c(1000, 1500, 1600, 1900, 2000, 2500, 3000, 3400,12000), c("#FCBBA1","#FC9C7D","#FB7D5E","#F85E42","#EF3B2C","#D32020","#B71318","#950B12" ,"#67000D"  ))
  
  continuous_lgd <- Legend(
    col_fun = col_fun,
    title = title,
    at = c(1000, 12000),
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontsize = 12),
    grid_height = unit(0.02, "cm"),
    grid_width = unit(0.50, "cm"),
    legend_height = unit(1.5, "cm")
  )
  
  discrete_lgd <- Legend(
    labels = names(fallback_colors),
    legend_gp = gpar(fill = fallback_colors),
    title = NULL,
    labels_gp = gpar(fontsize = 12),
    grid_height = unit(0.02, "cm"),
    grid_width = unit(0.50, "cm"),
    direction = "vertical"
  )
  
  packLegend(continuous_lgd, discrete_lgd, direction = "vertical")
}

make_centrifugation_legend2 <- function(title, range = c(1940, 16000), colors = c("#FDD0A2", "#7F2704")) {
  col_fun <- colorRamp2(c(1940, 2500, 3000, 6000, 12000, 13000, 15000, 16000), c("#FDD0A2","#FDB272", "#FD9649" ,"#F67824",     "#E65A0B","#CA4201","#A03303","#7F2704"))
  
  continuous_lgd <- Legend(
    col_fun = col_fun,
    title = title,
    at = c(1940, 16000),
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontsize = 12),
    grid_height = unit(0.02, "cm"),
    grid_width = unit(0.50, "cm"),
    legend_height = unit(1.5, "cm")
  )
  
  discrete_lgd <- Legend(
    labels = names(fallback_colors),
    legend_gp = gpar(fill = fallback_colors),
    title = NULL,
    labels_gp = gpar(fontsize = 12),
    grid_height = unit(0.02, "cm"),
    grid_width = unit(0.50, "cm"),
    direction = "vertical"
  )
  
  packLegend(continuous_lgd, discrete_lgd, direction = "vertical")
}


heatmap_list <- list()

for (var in row_order) {
  
  
  if (is.null(palette_list[[var]])) next
  
  values <- as.character(metadata_matrix[var, ])
  values[is.na(values)] <- "NA"
  
  if (var %in% c("centrifugation_step_1", "centrifugation_step_2")) {
    values <- factor(values, levels = c("1000g", "1500g", "1600g", "1900g", "1940g", "2000g", "2500g", "3000g", 
                                        "3400g", "6000g", "12000g", "13000g", "15000g", "16000g",
                                        "Unspecified", "placeholder", "None"))
    values <- as.character(values)
  } 
  if (var == "read_length") {
    values <- factor(values, levels = c("1x50", "1x75", "2x75", "2x100", "2x150", "NA"))
    values <- as.character(values)
  }
  
  color_map <- palette_list[[var]]
  if (is.null(color_map) || is.function(color_map)) next
  
  unique_vals <- unique(values)
  color_map <- color_map[unique_vals]
  names(color_map) <- unique_vals
  
  # Ensure specials
  if (!"NA" %in% names(color_map)) color_map["NA"] <- "grey80"
  color_map["Unspecified"] <- "grey20"
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
  } else if (var %in% c(
    "plasma_tubes_short_name", "biomaterial", "nucleic_acid_type",
    "rna_extraction_kit_short_name", "dnase", "library_prep_kit_short_name",
    "library_selection", "cdna_library_type", "broad_protocol_category"
  )) {
    special_levels <- c("Unspecified", "placeholder", "NA", "None")
    main_levels <- setdiff(names(color_map), special_levels)
    color_map <- color_map[c(sort(main_levels), intersect(special_levels, names(color_map)))]
  }
  
  special_levels <- c("Unspecified", "placeholder", "NA", "None")
  ordered_vals <- c(
    setdiff(names(color_map), special_levels),
    intersect(special_levels, names(color_map))
  )
  color_map <- color_map[ordered_vals]
  
  var_pretty <- clean_names[[var]]
  ordered_pretty_names <- clean_names[row_order]
  
  mat <- matrix(values, nrow = 1,
                dimnames = list(factor(var_pretty, levels = ordered_pretty_names),
                                colnames(metadata_matrix)))
  
  bottom_anno <- NULL
  if (var == "read_length") {
    bottom_anno <- HeatmapAnnotation(
      spacer = anno_empty(border = FALSE, height = unit(0.15, "cm")) # this controls the whitespace between RL and BPC rows
    )
  } else if (var =="broad_protocol_category") {
    bottom_anno <- HeatmapAnnotation(
      spacer = anno_empty(border = FALSE, height = unit(0.05, "cm")) # this controls the whitespace between BPC and column names, to give space to brackets
    )
  }
  
  ht_global_opt(font = "Arial", ADD = TRUE)
  ht_global_opt(
    heatmap_column_names_gp = gpar(fontsize = 12),
    heatmap_row_names_gp = gpar(fontsize = 12),
    legend_title_gp = gpar(fontsize = 12),
    legend_labels_gp = gpar(fontsize = 12, lineheight = 0.8)
  )
  
  ht <- Heatmap(
    mat,
    name = var_pretty,
    col = color_map,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_names_rot = 45,
    column_names_side = "bottom",
    column_names_gp = gpar(fontsize = 12),
    row_names_side = "left",
    bottom_annotation = bottom_anno,
    height = unit(0.32, "cm"), #controls the height of the heatmap "squares", this one is for the legend
    #height = unit(0.70, "cm"), #controls the height of the heatmap "squares", and this one for the figure
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.rect(x = x, y = y, width = width, height = height,
                gp = gpar(fill = fill, col = "white", lwd = 0.5))
    },
    heatmap_legend_param = list(
      grid_height = unit(0.02, "cm"),
      grid_width = unit(0.50, "cm"), # controls the shape of the legend squares
      gap = unit(0.1, "cm"),
      row_gap = unit(0.01, "cm"),
      title_gp = gpar(fontsize = 12, fontface = "bold"),
      labels_gp = gpar(fontsize = 12))
  )
  
  heatmap_list[[var]] <- ht
}



names(heatmap_list) <- row_order


ht_list <- Reduce(`%v%`, heatmap_list)




ragg::agg_png(
  "figures/fig_1b_metadata_heatmap_for_legend.png",
  width = 6.9, height = 4.6*0.766667, units = "in", res = 600, scaling = 5/12
)

# ragg::agg_png(
#   "figures/fig_1b_metadata_heatmap_for_figure.png",
#   width = 3.7, height = 2.5, units = "in", res = 600, scaling = 5/12
# )

ht_opt(legend_gap = unit(0.4, "mm"), ADD = TRUE) #reduce space between different legends



# ht_drw <- draw(ht_list,
#                # heatmap_legend_side = "right",
#                # annotation_legend_side = "right",
#                show_heatmap_legend = FALSE,
#                show_annotation_legend = FALSE,
#                gap = unit(0.2, "mm"))  # reduce from the default


ht_drw <- draw(
  ht_list,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  gap = unit(0.2, "mm"),
  heatmap_legend_list = list(
    "Centrifugation, step 1" = make_centrifugation_legend1("Centrifugation, step 1"),
    "Centrifugation, step 2" = make_centrifugation_legend2("Centrifugation, step 2")
  )
)


ht_pos <- htPositionsOnDevice(ht_drw)

y_top_in <- ht_pos[ht_pos$heatmap == "Broad protocol\ncategory (BPC)", "y_min"] - unit(0.45, "in") #this addition or subtraction here controls the position of the brackets. I could not find a better way to do it, for legend
#y_top_in <- ht_pos[ht_pos$heatmap == "Broad protocol\ncategory (BPC)", "y_min"] - unit(0.3, "in") #for figure

y_top_np <- convertY(y_top_in, "npc", valueOnly = FALSE)

y_bracket <- y_top_np + unit(0.1, "mm")          
y_vert_hi <- y_bracket
y_vert_lo <- y_bracket - unit(1, "mm")           

x_min_in <- ht_pos[1, "x_min"]
x_max_in <- ht_pos[1, "x_max"]
x_min_np <- convertX(x_min_in, "npc", valueOnly = TRUE)
x_max_np <- convertX(x_max_in, "npc", valueOnly = TRUE)
dx_np    <- x_max_np - x_min_np

seekViewport("global")

n_cols <- length(core_order)

for (i in seq_len(nrow(bracket_df))) {
  x1 <- bracket_df$xmin_idx[i]
  x2 <- bracket_df$xmax_idx[i]
  if (is.na(x1) || is.na(x2)) next
  
  offset <- 0.005
  x_start_np <- (x_min_np + (x1 - 1) / n_cols * dx_np) + offset
  x_end_np   <- (x_min_np + x2       / n_cols * dx_np) - offset
  
  x_start_u  <- unit(x_start_np, "npc") 
  x_end_u    <- unit(x_end_np,   "npc") 
  
  h_bar <- linesGrob(
    x = unit.c(x_start_u, x_end_u),
    y = unit.c(y_bracket, y_bracket),
    gp = gpar(col = "black", lwd = 1)
  )
  
  v_bars <- gList(
    linesGrob(x = unit.c(x_start_u, x_start_u),
              y = unit.c(y_vert_lo, y_vert_hi),
              gp = gpar(col = "black", lwd = 1)),
    linesGrob(x = unit.c(x_end_u, x_end_u),
              y = unit.c(y_vert_lo, y_vert_hi),
              gp = gpar(col = "black", lwd = 1))
  )
  
  grid.draw(grobTree(h_bar, v_bars))
}



dev.off()


svglite("figures/fig_1b_metadata_heatmap_for_legend.svg", width = 6.9, height = 4.6*0.766667, scaling = 5/12)
# svglite("figures/fig_1b_metadata_heatmap_for_figure.svg",
#          width = 3.7, height = 2.5, scaling = 5/12)
#showtext::showtext_begin()

# ht_drw <- draw(ht_list,
#                # heatmap_legend_side = "right",
#                # annotation_legend_side = "right",
#                show_heatmap_legend = FALSE,
#                show_annotation_legend = FALSE,
#                gap = unit(0.2, "mm"))  # reduce from the default

ht_drw <- draw(
  ht_list,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  gap = unit(0.2, "mm"),
  heatmap_legend_list = list(
    "Centrifugation, step 1" = make_centrifugation_legend1("Centrifugation, step 1"),
    "Centrifugation, step 2" = make_centrifugation_legend2("Centrifugation, step 2")
  )
)

ht_pos <- htPositionsOnDevice(ht_drw)

y_top_in <- ht_pos[ht_pos$heatmap == "Broad protocol\ncategory (BPC)", "y_min"] - unit(0.5, "in") #this addition or subtraction here controls the position of the brackets. I could not find a better way to do it
#y_top_in <- ht_pos[ht_pos$heatmap == "Broad protocol\ncategory (BPC)", "y_min"] - unit(0.3, "in") #for figure

y_top_np <- convertY(y_top_in, "npc", valueOnly = FALSE)

y_bracket <- y_top_np + unit(0.1, "mm")          
y_vert_hi <- y_bracket
y_vert_lo <- y_bracket - unit(1, "mm")           

x_min_in <- ht_pos[1, "x_min"]
x_max_in <- ht_pos[1, "x_max"]
x_min_np <- convertX(x_min_in, "npc", valueOnly = TRUE)
x_max_np <- convertX(x_max_in, "npc", valueOnly = TRUE)
dx_np    <- x_max_np - x_min_np

seekViewport("global")

n_cols <- length(core_order)

for (i in seq_len(nrow(bracket_df))) {
  x1 <- bracket_df$xmin_idx[i]
  x2 <- bracket_df$xmax_idx[i]
  if (is.na(x1) || is.na(x2)) next
  
  offset <- 0.005
  x_start_np <- (x_min_np + (x1 - 1) / n_cols * dx_np) + offset
  x_end_np   <- (x_min_np + x2       / n_cols * dx_np) - offset
  
  x_start_u  <- unit(x_start_np, "npc") 
  x_end_u    <- unit(x_end_np,   "npc") 
  
  h_bar <- linesGrob(
    x = unit.c(x_start_u, x_end_u),
    y = unit.c(y_bracket, y_bracket),
    gp = gpar(col = "black", lwd = 1)
  )
  
  v_bars <- gList(
    linesGrob(x = unit.c(x_start_u, x_start_u),
              y = unit.c(y_vert_lo, y_vert_hi),
              gp = gpar(col = "black", lwd = 1)),
    linesGrob(x = unit.c(x_end_u, x_end_u),
              y = unit.c(y_vert_lo, y_vert_hi),
              gp = gpar(col = "black", lwd = 1))
  )
  
  grid.draw(grobTree(h_bar, v_bars))
}



dev.off()

# variables to compare, comment out lines to play with combinations
preanalytical_vars <- c(
  "plasma_tubes_short_name",
  #"read_length",
  "centrifugation_step_1",
  "centrifugation_step_2",
  #"biomaterial",
  #"nucleic_acid_type",
  "rna_extraction_kit_short_name",
  "dnase",
  "library_prep_kit_short_name"#,
  #"library_selection",
  #"cdna_library_type"
)

protocol_combinations <- data_heatmap %>%
  select(dataset_batch, all_of(preanalytical_vars)) %>%
  mutate(across(everything(), ~ifelse(is.na(.), "NA", as.character(.)))) %>%
  distinct()  

shared_protocols <- protocol_combinations %>%
  group_by(across(all_of(preanalytical_vars))) %>%
  summarise(datasets = list(dataset_batch), n = n(), .groups = "drop") %>%
  arrange(desc(n))

print(shared_protocols)
shared_protocols$datasets
