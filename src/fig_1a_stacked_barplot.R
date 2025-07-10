library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(forcats)
library(jsonlite)
library(grid)
library(scales)
library(svglite)
#library(showtext)
#font_add(family="Arial", regular = "/usr/share/fonts/truetype/msttcorefonts/Arial.ttf")
#showtext_opts(dpi = 600)  # MUST come before showtext_auto()
#showtext_auto()
#theme_set(theme_classic(base_family = "Arial"))
library("extrafont")
loadfonts()

setwd("~/fl-cfRNAmeta/")

data <- read.delim("tables/cfRNA-meta_per_sample_metadata.tsv", header = TRUE, sep = "\t", fill = TRUE)

#Assign healthy phenotype to Flomics_1 and flomics_2
#data$phenotype[data$dataset_batch %in% c("flomics_1", "flomics_2", "giraldez_standard", "giraldez_phospho-rna-seq")] <- "healthy"

data_barplot <- data
data_barplot$simple_phenotype <- NA



data_barplot$simple_phenotype[data_barplot$phenotype == "healthy"] <- "healthy"
data_barplot$simple_phenotype[is.na(data_barplot$phenotype)] <- "missing"
data_barplot$simple_phenotype[data_barplot$phenotype == ""] <- "missing"
data_barplot$simple_phenotype[grep("Acute Myeloid Leukemia",data_barplot$phenotype)] <- "cancer"
data_barplot$simple_phenotype[grep("Alzheimers disease",data_barplot$phenotype)] <- "non-cancer disease"
data_barplot$simple_phenotype[grep("Chronic hepatitis B",data_barplot$phenotype)] <- "non-cancer disease"
data_barplot$simple_phenotype[grep("Chronic kidney failure EPO-treated",data_barplot$phenotype)] <- "non-cancer disease"
data_barplot$simple_phenotype[grep("Cirrhosis",data_barplot$phenotype)] <- "non-cancer disease"
data_barplot$simple_phenotype[grep("Colorectal cancer",data_barplot$phenotype)] <- "cancer"
data_barplot$simple_phenotype[grep("Diffuse large B-cell lymphoma",data_barplot$phenotype)] <- "cancer"
data_barplot$simple_phenotype[grep("Diverticulitis",data_barplot$phenotype)] <- "non-cancer disease"
data_barplot$simple_phenotype[grep("Esophagus cancer",data_barplot$phenotype)] <- "cancer"
data_barplot$simple_phenotype[grep("G-CSF-treated healthy donors",data_barplot$phenotype)] <- "healthy"
data_barplot$simple_phenotype[grep("Healthy",data_barplot$phenotype)] <- "healthy"
data_barplot$simple_phenotype[grep("Healthy pregnant woman",data_barplot$phenotype)] <- "healthy"
data_barplot$simple_phenotype[grep("Healthy pregnant woman who delivered preterm",data_barplot$phenotype)] <- "non-cancer disease"
data_barplot$simple_phenotype[grep("Liver cancer",data_barplot$phenotype)] <- "cancer"
data_barplot$simple_phenotype[grep("Lung cancer",data_barplot$phenotype)] <- "cancer"
data_barplot$simple_phenotype[grep("Multiple myeloma",data_barplot$phenotype)] <- "cancer"
data_barplot$simple_phenotype[grep("Nonalcoholic fatty liver disease",data_barplot$phenotype)] <- "non-cancer disease"
data_barplot$simple_phenotype[grep("Nonalcoholic steatohepatitis",data_barplot$phenotype)] <- "non-cancer disease"
data_barplot$simple_phenotype[grep("Pancreatic cancer",data_barplot$phenotype)] <- "cancer"
data_barplot$simple_phenotype[grep("Pre-cancerous condition: cirrhosis",data_barplot$phenotype)] <- "non-cancer disease"
data_barplot$simple_phenotype[grep("Pre-cancerous condition: MGUS",data_barplot$phenotype)] <- "non-cancer disease"
data_barplot$simple_phenotype[grep("Pre-eclampsia",data_barplot$phenotype)] <- "non-cancer disease"
data_barplot$simple_phenotype[grep("Primary mediastinal B-cell lymphoma",data_barplot$phenotype)] <- "cancer"
data_barplot$simple_phenotype[grep("Stomach cancer",data_barplot$phenotype)] <- "cancer"



mappings <- fromJSON("src/dataset_mappings.json")

clean_dataset_names <- unlist(mappings$datasetsLabels)
core_order <- unlist(mappings$datasetVisualOrder)

data_barplot$dataset_batch_clean <- recode(data_barplot$dataset_batch, !!!clean_dataset_names)

# all names alphabetically except "rozowsky" and "wei" last
ordered_names <- c(clean_dataset_names[core_order])

data_barplot$dataset_batch_clean <- factor(data_barplot$dataset_batch_clean, levels = ordered_names)

plot_data <- data_barplot %>%
  group_by(dataset_batch_clean, simple_phenotype) %>%
  summarise(count = n(), .groups = 'drop')

custom_colors <- brewer.pal(4, "Dark2")

# Option 1: simplified phenotype
ggplot(plot_data, aes(x = dataset_batch_clean, y = count, fill = simple_phenotype)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = brewer.pal(4, "Dark2"), name = "Simplified phenotype") +
  labs(x = "Dataset", y = "Number of samples") +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 18, color = "black"),
    axis.title.x = element_text(size = 20, face = "bold"), 
    axis.title.y = element_text(size = 20, face = "bold"),  
    legend.title = element_text(size = 18, face = "bold"),  
    legend.text = element_text(size = 16, face = "bold"),  
    panel.grid.major = element_line(size = 0.8),
    panel.grid.minor = element_blank(),
    plot.margin = margin(20, 20, 20, 20),
    plot.background = element_rect(
      fill = "white",
      colour = "white"
    )
  )

#ggsave("~/figures/fig_1a_sample_status_per_dataset.png", width = 15, height = 8, dpi = 600, units = "in")

# Option 2: detail into cancer type

# Filter original cancer cases only
cancer_detail_plot_data <- data_barplot %>%
  filter(simple_phenotype == "cancer") %>%
  group_by(dataset_batch_clean, phenotype) %>%
  summarise(count = n(), .groups = 'drop')

cancer_detail_plot_data$dataset_batch_clean <- factor(cancer_detail_plot_data$dataset_batch_clean,
                                                      levels = levels(data_barplot$dataset_batch_clean))

n_cancer_types <- length(unique(cancer_detail_plot_data$phenotype))
cancer_colors <- colorRampPalette(brewer.pal(10, "Set3"))(n_cancer_types)

ggplot(cancer_detail_plot_data, aes(x = dataset_batch_clean, y = count, fill = phenotype)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cancer_colors, name = "Cancer subtype") +
  labs(x = "Dataset", y = "Number of samples") +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 18, color = "black"),
    axis.title.x = element_text(size = 20, face = "bold"), 
    axis.title.y = element_text(size = 20, face = "bold"),  
    legend.title = element_text(size = 18, face = "bold"),  
    legend.text = element_text(size = 14),  
    panel.grid.major = element_line(size = 0.8),
    panel.grid.minor = element_blank(),
    plot.margin = margin(20, 20, 20, 20),
    plot.background = element_rect(fill = "white", colour = "white")
  )

#ggsave("~/figures/fig_1a_cancer_subtypes_per_dataset.png", width = 15, height = 8, dpi = 600, units = "in", bg = "white")

# Build third phenotype column: merge simple + detailed
data_barplot$phenotype_merged <- data_barplot$simple_phenotype
data_barplot$phenotype_merged[data_barplot$simple_phenotype == "cancer"] <- 
  data_barplot$phenotype[data_barplot$simple_phenotype == "cancer"]
data_barplot$phenotype_merged[data_barplot$phenotype_merged == "missing"] <- "Unspecified"
phenotype_merged_plot_data <- data_barplot %>%
  filter(!is.na(phenotype_merged)) %>%
  group_by(dataset_batch_clean, phenotype_merged) %>%
  summarise(count = n(), .groups = "drop")

merged_order <- c("healthy", "non-cancer disease", "Unspecified", sort(unique(data_barplot$phenotype[data_barplot$simple_phenotype == "cancer"])))
phenotype_merged_plot_data$phenotype_merged <- factor(phenotype_merged_plot_data$phenotype_merged, levels = merged_order)

n_cancer <- length(unique(data_barplot$phenotype[data_barplot$simple_phenotype == "cancer"]))
cancer_colors <- colorRampPalette(brewer.pal(10, "Set3"))(n_cancer)
names(cancer_colors) <- sort(unique(data_barplot$phenotype[data_barplot$simple_phenotype == "cancer"]))

merged_colors <- c(
  "Healthy" = "#1f78b4",
  "Non-cancer disease" = "#33a02c",
  "Unspecified" = "grey20",
  cancer_colors
)
bracket_df <- data.frame(
  xmin = c("Block (2x75bp)", "Giráldez (phospho)", "Ibarra (buffy coat)","Reggiardo (BioIVT)" , "Moufarrej (Site 1)", "Roskams (pilot)"),
  xmax = c("Block (2x150bp)", "Giráldez (standard)", "Ibarra (serum)","Reggiardo (DLS)" , "Moufarrej (Site 2)", "Roskams (validation)"),
  label = c("Block", "Giráldez", "Ibarra","Reggiardo" , "Moufarrej", "Roskams-Hieter")
)


phenotype_merged_plot_data$phenotype_merged <- fct_recode(
  phenotype_merged_plot_data$phenotype_merged,
  "Healthy" = "healthy",
  "Non-cancer disease" = "non-cancer disease"
)

legend_order <- c(
  "Healthy", "Non-cancer disease", "Liver cancer", "Colorectal cancer", "Esophagus cancer",
  "Lung cancer", "Stomach cancer", "Diffuse large B-cell lymphoma", "Primary mediastinal B-cell lymphoma",
  "Multiple myeloma","Acute Myeloid Leukemia" , "Pancreatic cancer", "Unspecified"
)

phenotype_merged_plot_data$phenotype_merged <- factor(phenotype_merged_plot_data$phenotype_merged, levels = legend_order)

merged_colors <- merged_colors[legend_order]

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



p <- ggplot(phenotype_merged_plot_data, aes(x = dataset_batch_clean, y = count, fill = phenotype_merged)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = merged_colors, name = "Donor phenotype") +
  labs(x = "Dataset", y = "Number of samples") +
  theme_minimal(base_size = 20) +
  theme(
    text=element_text(family="Arial"),
    axis.text.x = element_text(angle = 45, vjust = 1.1, hjust = 1, size = 12, color = "black"),
    axis.title.x = element_blank(), 
    axis.title.y = element_text(size = 20), 
    legend.position = "right",
    legend.title = element_text(size = 18, face = "bold"),  
    legend.text = element_text(size = 14),  
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(size = 0.8), 
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 40),
    plot.background = element_rect(fill = "white", colour = "white") 
  ) + coord_cartesian(clip = "off")#, ylim = c(-20, NA))




p <- add_bottom_brackets(p, bracket_df, levels(phenotype_merged_plot_data$dataset_batch_clean), y_base = 0.03)


ggsave("figures/fig_1a_combined_simplified_and_cancer_detail.png",
       width = 17, height = 8, dpi = 600, units = "in", bg = "white",device = ragg::agg_png)
ggsave("figures/fig_1a_combined_simplified_and_cancer_detail.svg",
       width = 17, height = 8, dpi = 600, units = "in", bg = "white", device = "svg")


library(ggplot2)
library(cowplot)
library(grid)  

plot_with_legend <- ggplot(phenotype_merged_plot_data, aes(x = dataset_batch_clean, y = count, fill = phenotype_merged)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = merged_colors, name = "Donor phenotype") +
  labs(x = "Dataset", y = "Number of samples") +
  theme_minimal(base_size = 6) +
  theme(
    text = element_text(family = "Arial", size = 12),
    axis.text.x = element_text(angle = 45, vjust = 1.05, hjust = 1, size = 12, color = "black"),
    axis.title.x = element_blank(),
    legend.key.size = unit(0.5, "cm"),  # or smaller, e.g., 0.25
    axis.title.y = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(size = 0.3),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(0, 0, 0, 27),
    plot.background = element_rect(fill = "white", colour = "white")
  ) +
  coord_cartesian(clip = "off")

plot_with_legend <- add_bottom_brackets(plot_with_legend, bracket_df, levels(phenotype_merged_plot_data$dataset_batch_clean), y_base = 0.03)


legend <- cowplot::get_legend(plot_with_legend)

plot_no_legend <- plot_with_legend + theme(legend.position = "none")
ggsave("figures/fig_1a_combined_simplified_and_cancer_detail_no_legend.png",
       plot = plot_no_legend,
       width = 3.35, height = 1.7, dpi = 600, units = "in", bg = "white", device = ragg::agg_png, scaling= 5/12)
ggsave("figures/fig_1a_combined_simplified_and_cancer_detail_legend_only.png",
       plot = legend,
       width = 1.3, height = 1.7, dpi = 600, units = "in", bg = "white", device = ragg::agg_png, scaling= 5/12)
ggsave("figures/fig_1a_combined_simplified_and_cancer_detail_no_legend.svg",
       plot = plot_no_legend,
       width = 3.35, height = 1.7, dpi = 600, units = "in", bg = "white", device = "svg", scaling= 5/12)
ggsave("figures/fig_1a_combined_simplified_and_cancer_detail_legend_only.svg",
       plot = legend,
       width = 1.3, height = 1.7, dpi = 600, units = "in", bg = "white", device = "svg", scaling= 5/12)


final_plot <- cowplot::plot_grid(
  plot_no_legend,
  legend,
  rel_widths = c(3.35, 1.3), 
  rel_heights = c(1,1),
  nrow = 1
)

ggsave("figures/fig_1a_combined_simplified_and_cancer_detail_2.png",
       plot = final_plot,
       width = 4.65, height = 1.7, dpi = 600, units = "in", bg = "white", device = ragg::agg_png, scaling = 5/12)

ggsave("figures/fig_1a_combined_simplified_and_cancer_detail_2.svg",
       plot = final_plot,
       width = 4.65, height = 1.7, dpi = 600, units = "in", bg = "white", device = "svg", scaling = 5/12)
