library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(forcats)
library(jsonlite)
library(grid)
library(scales)
library(showtext)
font_add("DejaVu Sans", regular = "DejaVuSans.ttf")
showtext_opts(dpi = 600)  # MUST come before showtext_auto()
showtext_auto()
theme_set(theme_classic(base_family = "DejaVu Sans"))



setwd("~/fl-cfRNAmeta/")

data <- read.table("tables/cfRNA-meta_per_sample_metadata.tsv", header = TRUE, sep = "\t", fill = TRUE)

#Assign healthy phenotype to Flomics_1 and flomics_2
data$phenotype[data$dataset_batch %in% c("flomics_1", "flomics_2", "giraldez_standard", "giraldez_phospho-rna-seq")] <- "healthy"

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
  "Unspecified" = "grey60",
  cancer_colors
)

phenotype_merged_plot_data$phenotype_merged <- fct_recode(
  phenotype_merged_plot_data$phenotype_merged,
  "Healthy" = "healthy",
  "Non-cancer disease" = "non-cancer disease"
)


ggplot(phenotype_merged_plot_data, aes(x = dataset_batch_clean, y = count, fill = phenotype_merged)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = merged_colors, name = "Phenotype / Cancer subtype") +
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

ggsave("~/figures/fig_1a_combined_simplified_and_cancer_detail.png",
       width = 15, height = 8, dpi = 600, units = "in", bg = "white",device = ragg::agg_png)
ggsave("~/figures/fig_1a_combined_simplified_and_cancer_detail.pdf",
       width = 15, height = 8, dpi = 600, units = "in", bg = "white", device = cairo_pdf)

