library(ggplot2)
library(dplyr)
library(ggpubr)
library(scales)
library(ggside)

################################################################################
# EXTERNAL DATASETS
################################################################################
# Read column names from text file
setwd("~/cfRNA-meta/")
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
                  "exonic_reads_minus_spike_ins"   )

data <- read.table("sampleinfo_snakeDA_2025-05-08.tsv", header = TRUE, sep = "\t", fileEncoding = "UTF-7")

biotype_data <- data[, 67:148]
biotype_data <- biotype_data %>% dplyr::select(-contains('_fc'))


biotype_data$total <- rowSums(biotype_data)
biotype_data$percent_of_reads_mapping_to_spike_ins <- biotype_data$spike_in / biotype_data$total
biotype_data$percent_of_reads_mapping_to_spike_ins <- biotype_data$percent_of_reads_mapping_to_spike_ins * 100

data$exonic_reads_minus_spike_ins <- data$exonic - data$spike_in
data$exonic_reads_minus_spike_ins <- (data$exonic_reads_minus_spike_ins / data$mapped_fragments)  * 100

#comparison <- data.frame (data$sample_id, data$sequencing_batch, data$exonic, biotype_data$total)
#comparison$difference <- comparison$data.exonic - comparison$biotype_data.total
#summary(comparison$difference)
# write.csv(comparison, file = "exonic_minus_biotype.csv")

# Keep specified columns
selected_columns <- NULL
selected_columns <- c( "sample_id", "sequencing_batch", "status", "spike_in_pct", "protein_coding_pct", "percentage_of_spliced_reads", column_names)
filtered_data <- data[ ,selected_columns]

filtered_data$percent_of_reads_mapping_to_spike_ins <- biotype_data$percent_of_reads_mapping_to_spike_ins

result <- filtered_data %>%
  filter(percent_of_reads_mapping_to_spike_ins > 5) %>% 
  group_by(sequencing_batch) %>%     
  summarise(count = n())    

print(result)
# table_filtered <- filtered_data %>%
#   filter(percent_of_reads_mapping_to_spike_ins <= 5)
table_filtered <- filtered_data

table_filtered$log_genes_80 <- log(table_filtered$genes_contributing_to_80._of_reads)

# Generate colors from the glasbey palette
num_datasets <- length(unique(table_filtered$sequencing_batch))
glasbey_colors <- pals::glasbey(num_datasets)
table_filtered$sequencing_batch <- factor(table_filtered$sequencing_batch, levels = c("block_1",
                                                                                      "block_2", 
                                                                                      "decruyenaere",
                                                                                      "zhu",
                                                                                      "chen",
                                                                                      "ngo",
                                                                                      "roskams_1",
                                                                                      "roskams_2",
                                                                                      "moufarrej",
                                                                                      "sun_1",
                                                                                      "sun_2",
                                                                                      "tao",
                                                                                      "toden",
                                                                                      "ibarra",
                                                                                      "chalasani",
                                                                                      "rozowsky",
                                                                                      "taowei",
                                                                                      "giraldez",
                                                                                      "reggiardo",
                                                                                      "wang"))

datasetsPalette=c( "block_1" = "#b3b3b3",
                   "block_2" = "#7d7a7a",
                   "decruyenaere" =  "#009E73",
                   "zhu" = "#ffd633",
                   "chen" = "#997a00",
                   "ngo" =  "#fa8072",
                   "roskams_1" = "#944dff",
                   "roskams_2" = "#5b2e9e",
                   "moufarrej" = "#CC79A7",
                   "sun_1" = "#D55E00",
                   "sun_2" = "#8a3d00", 
                   "tao" ="#0072B2",
                   "toden" = "#800099",
                   "ibarra" = "#800000",
                   "chalasani" = "#800040",
                   "rozowsky" = "#006600",
                   "taowei"="#B32400",
                   "giraldez" = "#B1CC71",
                   "reggiardo" = "#F1085C",
                   "wang" = "#FE8F42") 

datasetsLabels=c("block_1" = "Block_1",
                 "block_2" = "Block_2",
                 "chen" = "Chen",
                 "decruyenaere" = "Decruyenaere",
                 "ngo" = "Ngo",
                 "roskams_1" = "Roskams-Hieter_1",
                 "roskams_2" = "Roskams-Hieter_2",
                 "moufarrej" = "Moufarrej",
                 "toden" = "Toden",
                 "ibarra" = "Ibarra",
                 "chalasani" = "Chalasani",
                 "rozowsky"="ENCODE\n(bulk tissue RNA-Seq)",
                 "sun_1" = "Sun_1",
                 "sun_2" = "Sun_2",
                 "tao" = "Tao",
                 "zhu" ="Zhu",
                 "taowei" = "Wei (cfDNA)",
                 "giraldez" = "Giráldez",
                 "reggiardo" = "Reggiardo",
                 "wang" = "Wang")

column_names <- c(column_names, "percent_of_reads_mapping_to_spike_ins", "log_genes_80")

#write.table(table_filtered, file="Qc_table_filtered.tsv", row.names = FALSE)
ggplot_objects <- lapply(column_names, function(col_name) {
  ggplot(table_filtered, aes(x = sequencing_batch, y = .data[[col_name]], fill = sequencing_batch)) +
    geom_boxplot(alpha = 0.3, color = "black", position = position_dodge(width = 0.75), outlier.shape = NA ) +
    geom_point(aes(y = .data[[col_name]], color = sequencing_batch), 
               position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.8), 
               shape = 16, size = 2) +
    labs(title = col_name,
         x = "Dataset", y = col_name) +
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
  #guides(fill = guide_legend(override.aes = list(color = "black"), ncol = 1)) 
  
})

setwd("~/cfRNA-meta/full_comparison_2025_05_08/")

# Save individual  boxplots 
for (i in 1:length(ggplot_objects)) {
  col_name <- column_names[i]
  output_file <- paste0(gsub(" ", "_", tolower(col_name)), "_external_datasets_boxplot_with_points.png")
  ggsave(output_file, ggplot_objects[[i]], width = 9, height = 6, dpi = 600)
  output_file <- paste0(gsub(" ", "_", tolower(col_name)), "_external_datasets_boxplot_with_points.pdf")
  ggsave(output_file, ggplot_objects[[i]], width = 9, height = 6, dpi = 600)
}


# p <- ggplot(data = table_filtered, aes(x = Exonic_percentage, y = log(Genes_contributing_to_80._of_reads), color = dataset)) +
#   geom_point() +
#   geom_density_2d() + # comment to remove the density lines
#   labs(x = "exonic percentage",
#        y = "log(# of genes contributing to 80% of reads)",
#        color = "") +
#   theme_minimal() +
#   scale_color_manual(values = datasetsPalette, labels = datasetsLabels)+
#   scale_fill_manual(values = datasetsPalette) 


########### Diversity scatterplot

library(scales)
p <- ggplot(data = table_filtered, aes(x = exonic_reads_minus_spike_ins, y = genes_contributing_to_80._of_reads, color = sequencing_batch)) +
  geom_point(size = 2, alpha = 0.7) + 
  geom_density_2d(aes(color = sequencing_batch), alpha = 0.5, size = 0.8) +
  labs(title = "",
       x = "Exonic percentage minus spike-ins",
       y = "# of Genes Contributing to 80% of Reads",
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

ggsave("diversity_scatterplot.png", plot = ps, width = 12, height = 8, dpi = 300)
ggsave("diversity_scatterplot.pdf", plot = ps, width = 12, height = 8, dpi = 300)

################## Shannon entropy

shannon <- read.table("~/cfRNA-meta/exp_mat/shannon_entropy.csv", header = TRUE, sep = ",")
colnames(shannon) <- c("sequencing_batch", "sample_id", "ShannonEntropy")

table_filtered_combined <- left_join(table_filtered, shannon, by = "sample_id")

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
  geom_boxplot(alpha = 0.3, color = "black", position = position_dodge(width = 0.75), outlier.shape = NA ) +
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
  scale_color_manual(values = datasetsPalette, labels = datasetsLabels)

ggsave("shannon_entropy_boxplot.png", plot = p, width = 9, height = 6, dpi = 600)
ggsave("shannon_entropy_boxplot.pdf", plot = p, width = 9, height = 6, dpi = 600)


############### Fragments mapped to the expected strand


table_filtered$Fragments_mapping_to_expected_strand_pct <- 100 - as.numeric(table_filtered$reads_mapping_sense_percentage)

p <- ggplot(table_filtered, aes(x = sequencing_batch, y = Fragments_mapping_to_expected_strand_pct, fill = sequencing_batch)) +
  geom_boxplot(alpha = 0.6, color = "black", position = position_dodge(width = 0.75), outlier.shape = NA) +
  geom_point(aes(y = Fragments_mapping_to_expected_strand_pct, color = sequencing_batch), 
             position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.8), 
             shape = 21, size = 1.5, stroke = 0.2, alpha = 0.6) +
  geom_hline(yintercept=100, linetype='dashed', col = 'lightgrey')+
  geom_hline(yintercept=50, linetype='dashed', col = 'lightgrey')+
  labs(title = "",
       x = "Dataset", y = "% fragments mapping to correct gene orientation") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "none") +
  scale_x_discrete(labels = datasetsLabels) +
  scale_fill_manual(values = datasetsPalette, labels = datasetsLabels) +
  scale_color_manual(values = datasetsPalette, labels = datasetsLabels)


ggsave("fragments_mapped_expected_strand.png", p, width = 9, height = 6, dpi = 600)
ggsave("fragments_mapped_expected_strand.pdf", p, width = 9, height = 6, dpi = 600)


################### Fragment number


table_filtered <- table_filtered %>%
  mutate(fragment_number = if_else(
    `sequencing_batch` == "giraldez",
    read_number,           # Single-end: keep as is
    read_number / 2        # Paired-end: divide by 2
  ))

p <- ggplot(table_filtered, aes(x = sequencing_batch, y = fragment_number, fill = sequencing_batch)) +
  geom_boxplot(alpha = 0.6, color = "black", position = position_dodge(width = 0.75), outlier.shape = NA) +
  geom_point(aes(y = fragment_number, color = sequencing_batch), 
             position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.8), 
             shape = 21, size = 1.5, stroke = 0.2, alpha = 0.6) +
  labs(title = "",
       x = "Dataset", y = "Fragment number") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "none") +
  scale_x_discrete(labels = datasetsLabels) +
  scale_fill_manual(values = datasetsPalette, labels = datasetsLabels) +
  scale_color_manual(values = datasetsPalette, labels = datasetsLabels)

ggsave("fragment_number.png", p, width = 9, height = 6)
ggsave("fragment_number.pdf", p, width = 9, height = 6)


################################################################################
# Combine with Flomics liquidx for META-ANALYSIS
################################################################################

# Read column names from text file
setwd("~/cfRNA-meta/")
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
                  "exonic_reads_minus_spike_ins"   )

data <- read.table("sampleinfo_snakeDA_2025-05-08.tsv", header = TRUE, sep = "\t", fileEncoding = "UTF-7")
data_2nd_gen <- ""
data_liquidx <- read.table("sampleinfo_snakeDA_liquidx_2025_05_08.tsv", header = TRUE, sep = "\t", fileEncoding = "UTF-8")

# Remove inserm CRC
data_liquidx <- data_liquidx[data_liquidx$collection_subcenter_health_care_unit != "INSERM-CRC23",]
#Remove non-healthy samples (for meta-analysis only)
#data_liquidx <- data_liquidx[data_liquidx$status == "healthy"]
# Remove samples not passing QC
#data_liquidx <- data_liquidx[data_liquidx$percentage_of_spliced_reads > 20,]

### Merge the two dataframes
merged_df <- merge(data, data_liquidx, all = TRUE)


biotype_data <- merged_df[, 66:147]
biotype_data <- biotype_data %>% dplyr::select(-contains('_fc'))

biotype_data$total <- rowSums(biotype_data)
biotype_data$percent_of_reads_mapping_to_spike_ins <- biotype_data$spike_in / biotype_data$total
biotype_data$percent_of_reads_mapping_to_spike_ins <- biotype_data$percent_of_reads_mapping_to_spike_ins * 100

merged_df$exonic_reads_minus_spike_ins <- merged_df$exonic - merged_df$spike_in
merged_df$exonic_reads_minus_spike_ins <- (merged_df$exonic_reads_minus_spike_ins / merged_df$mapped_fragments)  * 100

#comparison <- data.frame (data$sample_id, data$sequencing_batch, data$exonic, biotype_data$total)
#comparison$difference <- comparison$data.exonic - comparison$biotype_data.total
#summary(comparison$difference)
# write.csv(comparison, file = "exonic_minus_biotype.csv")

# Keep specified columns
selected_columns <- NULL
selected_columns <- c( "sample_id", "sequencing_batch", "status", "spike_in_pct", "protein_coding_pct", "percentage_of_spliced_reads", column_names)
filtered_data <- merged_df[ ,selected_columns]

filtered_data$percent_of_reads_mapping_to_spike_ins <- biotype_data$percent_of_reads_mapping_to_spike_ins

result <- filtered_data %>%
  filter(percent_of_reads_mapping_to_spike_ins > 5) %>% 
  group_by(sequencing_batch) %>%     
  summarise(count = n())    

print(result)
# table_filtered <- filtered_data %>%
#   filter(percent_of_reads_mapping_to_spike_ins <= 5)
table_filtered <- filtered_data

table_filtered$log_genes_80 <- log(table_filtered$genes_contributing_to_80._of_reads)

# transform all "FL-" sequencing batches to "liquidx"
table_filtered$sequencing_batch[grepl("^FL", table_filtered$sequencing_batch)] <- "liquidx"

num_datasets <- length(unique(table_filtered$sequencing_batch))
glasbey_colors <- pals::glasbey(num_datasets)
table_filtered$sequencing_batch <- factor(table_filtered$sequencing_batch, levels = c("liquidx",
                                                                                      "block_1",
                                                                                      "block_2", 
                                                                                      "decruyenaere",
                                                                                      "zhu",
                                                                                      "chen",
                                                                                      "ngo",
                                                                                      "roskams_1",
                                                                                      "roskams_2",
                                                                                      "moufarrej",
                                                                                      "sun_1",
                                                                                      "sun_2",
                                                                                      "tao",
                                                                                      "toden",
                                                                                      "ibarra",
                                                                                      "chalasani",
                                                                                      "rozowsky",
                                                                                      "taowei",
                                                                                      "giraldez",
                                                                                      "reggiardo",
                                                                                      "wang"))

datasetsPalette=c( "liquidx" = "#144d6b", 
                   "block_1" = "#b3b3b3",
                   "block_2" = "#7d7a7a",
                   "decruyenaere" =  "#009E73",
                   "zhu" = "#ffd633",
                   "chen" = "#997a00",
                   "ngo" =  "#fa8072",
                   "roskams_1" = "#944dff",
                   "roskams_2" = "#5b2e9e",
                   "moufarrej" = "#CC79A7",
                   "sun_1" = "#D55E00",
                   "sun_2" = "#8a3d00", 
                   "tao" ="#0072B2",
                   "toden" = "#800099",
                   "ibarra" = "#800000",
                   "chalasani" = "#800040",
                   "rozowsky" = "#006600",
                   "taowei"="#B32400",
                   "giraldez" = "#B1CC71",
                   "reggiardo" = "#F1085C",
                   "wang" = "#FE8F42") 

datasetsLabels=c("liquidx" = "Flomics",
                 "block_1" = "Block_1",
                 "block_2" = "Block_2",
                 "chen" = "Chen",
                 "decruyenaere" = "Decruyenaere",
                 "ngo" = "Ngo",
                 "roskams_1" = "Roskams-Hieter_1",
                 "roskams_2" = "Roskams-Hieter_2",
                 "moufarrej" = "Moufarrej",
                 "toden" = "Toden",
                 "ibarra" = "Ibarra",
                 "chalasani" = "Chalasani",
                 "rozowsky"="ENCODE\n(bulk tissue RNA-Seq)",
                 "sun_1" = "Sun_1",
                 "sun_2" = "Sun_2",
                 "tao" = "Tao",
                 "zhu" ="Zhu",
                 "taowei" = "Wei (cfDNA)",
                 "giraldez" = "Giráldez",
                 "reggiardo" = "Reggiardo",
                 "wang" = "Wang")

column_names <- c(column_names, "percent_of_reads_mapping_to_spike_ins", "log_genes_80")

#write.table(table_filtered, file="Qc_table_filtered.tsv", row.names = FALSE)
ggplot_objects <- lapply(column_names, function(col_name) {
  ggplot(table_filtered, aes(x = sequencing_batch, y = .data[[col_name]], fill = sequencing_batch)) +
    geom_boxplot(alpha = 0.3, color = "black", position = position_dodge(width = 0.75), outlier.shape = NA ) +
    geom_point(aes(y = .data[[col_name]], color = sequencing_batch), 
               position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.8), 
               shape = 16, size = 2) +
    labs(title = col_name,
         x = "Dataset", y = col_name) +
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
  #guides(fill = guide_legend(override.aes = list(color = "black"), ncol = 1))
  
})



setwd("~/cfRNA-meta/full_comparison_2025_05_08/")

for (i in 1:length(ggplot_objects)) {
  col_name <- column_names[i]
  output_file <- paste0(gsub(" ", "_", tolower(col_name)), "_external_datasets_boxplot_with_points.png")
  ggsave(output_file, ggplot_objects[[i]], width = 9, height = 6, dpi = 600)
  output_file <- paste0(gsub(" ", "_", tolower(col_name)), "_external_datasets_boxplot_with_points.pdf")
  ggsave(output_file, ggplot_objects[[i]], width = 9, height = 6, dpi = 600)
}


# p <- ggplot(data = table_filtered, aes(x = Exonic_percentage, y = log(Genes_contributing_to_80._of_reads), color = dataset)) +
#   geom_point() +
#   geom_density_2d() + # comment to remove the density lines
#   labs(x = "exonic percentage",
#        y = "log(# of genes contributing to 80% of reads)",
#        color = "") +
#   theme_minimal() +
#   scale_color_manual(values = datasetsPalette, labels = datasetsLabels)+
#   scale_fill_manual(values = datasetsPalette) 


########### Diversity scatterplot

library(scales)
p <- ggplot(data = table_filtered, aes(x = exonic_reads_minus_spike_ins, y = genes_contributing_to_80._of_reads, color = sequencing_batch)) +
  geom_point(size = 2, alpha = 0.7) + 
  geom_density_2d(aes(color = sequencing_batch), alpha = 0.5, size = 0.8) +
  labs(title = "",
       x = "Exonic percentage minus spike-ins",
       y = "# of Genes Contributing to 80% of Reads",
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

ggsave("diversity_scatterplot.png", plot = ps, width = 12, height = 8, dpi = 300)
ggsave("diversity_scatterplot.pdf", plot = ps, width = 12, height = 8, dpi = 300)

################## Shannon entropy

shannon <- read.table("~/cfRNA-meta/exp_mat/shannon_entropy.csv", header = TRUE, sep = ",")
colnames(shannon) <- c("sequencing_batch", "sample_id", "ShannonEntropy")

table_filtered_combined <- left_join(table_filtered, shannon, by = "sample_id")

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
  geom_boxplot(alpha = 0.3, color = "black", position = position_dodge(width = 0.75), outlier.shape = NA ) +
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
  scale_color_manual(values = datasetsPalette, labels = datasetsLabels)

ggsave("shannon_entropy_boxplot.png", plot = p, width = 9, height = 6, dpi = 600)
ggsave("shannon_entropy_boxplot.pdf", plot = p, width = 9, height = 6, dpi = 600)


############### Fragments mapped to the expected strand


table_filtered$Fragments_mapping_to_expected_strand_pct <- 100 - as.numeric(table_filtered$reads_mapping_sense_percentage)

p <- ggplot(table_filtered, aes(x = sequencing_batch, y = Fragments_mapping_to_expected_strand_pct, fill = sequencing_batch)) +
  geom_boxplot(alpha = 0.6, color = "black", position = position_dodge(width = 0.75), outlier.shape = NA) +
  geom_point(aes(y = Fragments_mapping_to_expected_strand_pct, color = sequencing_batch), 
             position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.8), 
             shape = 21, size = 1.5, stroke = 0.2, alpha = 0.6) +
  geom_hline(yintercept=100, linetype='dashed', col = 'lightgrey')+
  geom_hline(yintercept=50, linetype='dashed', col = 'lightgrey')+
  labs(title = "",
       x = "Dataset", y = "% fragments mapping to correct gene orientation") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "none") +
  scale_x_discrete(labels = datasetsLabels) +
  scale_fill_manual(values = datasetsPalette, labels = datasetsLabels) +
  scale_color_manual(values = datasetsPalette, labels = datasetsLabels)


ggsave("fragments_mapped_expected_strand.png", p, width = 9, height = 6, dpi = 600)
ggsave("fragments_mapped_expected_strand.pdf", p, width = 9, height = 6, dpi = 600)


################### Fragment number


table_filtered <- table_filtered %>%
  mutate(fragment_number = if_else(
    `sequencing_batch` == "giraldez",
    read_number,           # Single-end: keep as is
    read_number / 2        # Paired-end: divide by 2
  ))

p <- ggplot(table_filtered, aes(x = sequencing_batch, y = fragment_number, fill = sequencing_batch)) +
  geom_boxplot(alpha = 0.6, color = "black", position = position_dodge(width = 0.75), outlier.shape = NA) +
  geom_point(aes(y = fragment_number, color = sequencing_batch), 
             position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.8), 
             shape = 21, size = 1.5, stroke = 0.2, alpha = 0.6) +
  labs(title = "",
       x = "Dataset", y = "Fragment number") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "none") +
  scale_x_discrete(labels = datasetsLabels) +
  scale_fill_manual(values = datasetsPalette, labels = datasetsLabels) +
  scale_color_manual(values = datasetsPalette, labels = datasetsLabels)

ggsave("fragment_number.png", p, width = 9, height = 6)
ggsave("fragment_number.pdf", p, width = 9, height = 6)



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


p <- ggplot(table_filtered_slide_deck, aes(x = dataset, y = percentage_of_spliced_reads, fill = dataset)) +
  geom_boxplot(alpha = 0.6, color = "black", position = position_dodge(width = 0.75), outlier.shape = NA) +
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
  scale_color_manual(values = datasetsPalette, labels = datasetsLabels) + # Match outline colors to fill
  guides(fill = guide_legend(override.aes = list(color = "black"), ncol = 1)) 

ggsave("percent_of_spliced_reads_slide_deck.png", p, width = 8, height = 6)
ggsave("percent_of_spliced_reads_slide_deck.pdf", p, width = 8, height = 6)


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
    log_lib_diversity = log2(genes_contributing_to_80._of_reads)
  ) %>%
  filter(!is.na(percentage_of_spliced_reads), !is.na(log_lib_diversity))



p <- ggplot(table_plot, aes(x = percentage_of_spliced_reads, y = log_lib_diversity, color = dataset)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_density_2d(alpha = 0.5, size = 0.8) +
  
  geom_xsideboxplot(
    aes(y = dataset, fill = dataset), orientation = "y",
    #width = 0.4,
    outlier.shape = 1,
    alpha = 0.3,
    #size = 0.4,
    side = "top",
    show.legend = FALSE
  ) +
  
  geom_ysideboxplot(
    aes(x = dataset, group = dataset, fill = dataset), orientation = "x",
    width = 0.4,
    outlier.shape = 1,
    alpha = 0.3,
    size = 0.4,
    side = "right",
    show.legend = FALSE
  ) +
  
  scale_color_manual(values = datasetsPalette, labels = datasetsLabels) +
  scale_fill_manual(values = datasetsPalette) +
  
  labs(
    x = "Percent of spliced reads",
    y = "Log2(Library diversity)",
    color = ""
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "left",
    ggside.panel.scale = 0.3, 
    ggside.axis.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    ggside.axis.text.x = element_text(angle = 45, hjust = 1),
    plot.background = element_rect(fill = "white", colour = "white")
  )

ggsave("library_vs_spliced_boxplots.png", plot = p, width = 17, height = 8, dpi = 600)


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


################################################################################
# Combined boxplots
################################################################################
# Load the "old" combined table with Flomics Generations
data <- read.table("~/cfRNA-meta/QC_table_flomics_gen_1-4.csv", header = TRUE, sep = "\t")
column_names <- readLines("~/cfRNA-meta/columns_old_table.txt")
column_names <- gsub("%", ".", column_names)

biotype_data <- data[, 50:115]
biotype_data <- biotype_data %>% dplyr::select(-contains('_fc'))

biotype_data$total <- rowSums(biotype_data)
biotype_data$percent_of_reads_mapping_to_spike_ins <- biotype_data$spike_in / biotype_data$total
biotype_data$percent_of_reads_mapping_to_spike_ins <- biotype_data$percent_of_reads_mapping_to_spike_ins * 100

data$exonic_reads_minus_spike_ins <- data$Exonic - data$spike_in
data$exonic_reads_minus_spike_ins <- (data$exonic_reads_minus_spike_ins / data$uniquely_mapped_reads)  * 100

# comparison <- data.frame (data$sample_id, data$sequencing_batch, data$exonic, biotype_data$total)
# comparison$difference <- comparison$data.exonic - comparison$biotype_data.total
# summary(comparison$difference)
# write.csv(comparison, file = "exonic_minus_biotype.csv")

# Keep specified columns
selected_columns <- NULL
selected_columns <- c( "Sample", "dataset", column_names)
filtered_data_old <- data[ ,selected_columns]

filtered_data_old$percent_of_reads_mapping_to_spike_ins <- biotype_data$percent_of_reads_mapping_to_spike_ins

result <- filtered_data_old %>%
  filter(percent_of_reads_mapping_to_spike_ins > 5) %>%  # Filter rows where variable1 is greater than 5
  group_by(dataset) %>%     # Group data by the 'dataset' column
  summarise(count = n())    # Count the number of rows in each group

print(result)
table_filtered_old <- filtered_data_old %>%
  filter(percent_of_reads_mapping_to_spike_ins <= 5)
table_filtered_old <- filtered_data_old

table_filtered_old$log_genes_80 <- log(table_filtered_old$number_of_genes_contributing_to_80._of_reads)

##############################################################################
# OLD CODE, DISREGARD FOR THE MOMENT (FOR COMPARATIVE META-ANALYSIS)
###############################################################################

# joint_table_final <- joint_table %>%
#   bind_rows(table_filtered_old)
# 
# # Add flomics_gen_4 to liquidx
# #joint_table_final[joint_table_final$dataset=="flomics_gen_4",] <- "liquidx"
# #joint_table_final$dataset <- gsub("flomics_gen_4", "liquidx", joint_table_final$dataset)
# #remove flomics_gen_1
# joint_table_final <- joint_table_final[joint_table_final$dataset != "flomics_gen_1", ]
# joint_table_final <- joint_table_final[joint_table_final$dataset != "flomics_gen_4", ]
# 
# 
# # Remove some datasets for SLIDE DECK
# joint_table_final <- joint_table_final[joint_table_final$dataset != "flomics_gen_1", ]
# joint_table_final <- joint_table_final[joint_table_final$dataset != "flomics_gen_2", ]
# joint_table_final <- joint_table_final[joint_table_final$dataset != "flomics_gen_3", ]
# joint_table_final <- joint_table_final[joint_table_final$dataset != "flomics_gen_4", ]
# joint_table_final <- joint_table_final[joint_table_final$dataset != "chalasani", ]
# joint_table_final <- joint_table_final[joint_table_final$dataset != "toden", ]
# joint_table_final <- joint_table_final[joint_table_final$dataset != "ibarra", ]
# joint_table_final <- joint_table_final[joint_table_final$dataset != "sun", ]
# joint_table_final <- joint_table_final[joint_table_final$dataset != "reggiardo", ]
# joint_table_final <- joint_table_final[joint_table_final$dataset != "wang", ]
# joint_table_final <- joint_table_final[joint_table_final$dataset != "giraldez", ]
# 
# # Generate colors from the glasbey palette
# num_datasets <- length(unique(joint_table_final$dataset))
# glasbey_colors <- pals::glasbey(num_datasets)
# 
# 
# 
# joint_table_final$dataset <- factor(joint_table_final$dataset, levels = c("flomics_gen_2","flomics_gen_3" ,"liquidx","block" ,"decruyenaere", "zhu", "chen",  "ngo", "roskams", "moufarrej","Sun","tao", "toden", "ibarra", "chalasani","rozowsky", "taowei"))
# joint_table_final$dataset <- factor(joint_table_final$dataset, levels = c( "liquidx","block" ,"decru", "zhu", "chen",  "ngo", "roskams", "moufarrej","tao","rozowsky", "taowei"))
# 
# datasetsPalette=c( "liquidx" = "#144d6b", "flomics_gen_2" = "#9AB9D6", "flomics_gen_3" = "#A5CAE5","block" = "#b3b3b3" , "decruyenaere" =  "#009E73" , "zhu" = "#ffd633", "chen" = "#997a00", "ngo" =  "#fa8072", "roskams" = "#944dff" , "moufarrej" = "#CC79A7", "Sun" = "#D55E00", "tao" ="#0072B2", "toden" = "#800099",  "ibarra" = "#800000", "chalasani" = "#800040","rozowsky" = "#006600", "taowei"="#B32400")
# datasetsPalette=c( "liquidx" = "#144d6b", "block" = "#b3b3b3" , "decru" =  "#009E73" , "zhu" = "#ffd633", "chen" = "#997a00", "ngo" =  "#fa8072", "roskams" = "#944dff" , "moufarrej" = "#CC79A7",  "tao" ="#0072B2","rozowsky" = "#006600", "taowei"="#B32400")
# 
# datasetsLabels=c("liquidx" = "Flomics 2024-Q3", "flomics_gen_2"="Flomics 2023-Q1", "flomics_gen_3"="Flomics 2023-Q4","block" = "Block 2022", "chen" = "Lu 2022", "decruyenaere" = "Decruyenaere 2023", "ngo" = "Ngo 2018", "roskams" = "Roskams 2022","moufarrej" = "Moufarrej 2022", "toden" = "Toden 2020", "ibarra" = "Ibarra 2020", "chalasani" = "Chalasani 2021","rozowsky"="ENCODE 2023\n(bulk tissue RNA-Seq)", "Sun" = "Sun 2024", "tao" = "Lu 2023", "zhu" ="Lu 2021", "taowei" = "Tao Wei (cfDNA)")
# datasetsLabels=c("liquidx" = "Flomics 2024-Q3", "block" = "Block 2022", "chen" = "Lu 2022", "decru" = "Decruyenaere 2023", "ngo" = "Ngo 2018", "roskams" = "Roskams 2022","moufarrej" = "Moufarrej 2022","rozowsky"="ENCODE 2023\n(bulk tissue RNA-Seq)", "tao" = "Lu 2023", "zhu" ="Lu 2021", "taowei" = "Tao Wei (cfDNA)")
# 
# #column_names <- c(column_names, "percent_of_reads_mapping_to_spike_ins", "log_genes_80")
# 
# 
# genomic_sun_samples <- c("SRR20968811", "SRR20968812", "SRR20968813", "SRR20968814")
# joint_table_final <- joint_table_final[!(joint_table_final$Sample %in% genomic_sun_samples),]
# 
# 
# joint_table_final$Fragments_mapping_to_expected_strand_pct <- 100 - as.numeric(joint_table_final$Reads_mapping_sense_percentage)
# 
# column_names <- c(column_names,"Fragments_mapping_to_expected_strand_pct")
# 
# #write.table(table_filtered, file="Qc_table_filtered.tsv", row.names = FALSE)
# library(ggplot2)
# 
# column_names <- readLines("~/cfRNA-meta/columns.txt")
# column_names <- gsub("%", ".", column_names)
# column_names <- c(column_names, "percent_of_reads_mapping_to_spike_ins"   , "log_genes_80"  )
# 
# ggplot_objects <- lapply(column_names, function(col_name) {
#   ggplot(joint_table_final, aes(x = dataset, y = .data[[col_name]], fill = dataset)) +
#     geom_boxplot(alpha = 0.6, color = "black", position = position_dodge(width = 0.75), outlier.shape = NA) +
#     geom_point(aes(y = .data[[col_name]], color = dataset), 
#                position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.8), 
#                shape = 21, size = 1.5, stroke = 0.2, alpha = 0.6) +
#     labs(title = col_name,
#          x = "Dataset", y = col_name) +
#     theme_classic() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
#           axis.title = element_text(size = 12, face = "bold"),
#           plot.title = element_text(size = 14, face = "bold"),
#           legend.position = "none", 
#           legend.title = element_blank(),
#           legend.key.size = unit(0.5, "cm"),
#           legend.spacing.y = unit(0.2, 'cm')) +
#     scale_x_discrete(labels = datasetsLabels) +
#     scale_fill_manual(values = datasetsPalette, labels = datasetsLabels) +
#     scale_color_manual(values = datasetsPalette, labels = datasetsLabels) + # Match outline colors to fill
#     guides(fill = guide_legend(override.aes = list(color = "black"), ncol = 1)) # Adjust legend for compactness
# })
# 
# p <- ggplot(joint_table_final, aes(x = dataset, y = Fragments_mapping_to_expected_strand_pct, fill = dataset)) +
#   geom_boxplot(alpha = 0.6, color = "black", position = position_dodge(width = 0.75), outlier.shape = NA) +
#   geom_point(aes(y = Fragments_mapping_to_expected_strand_pct, color = dataset), 
#              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.8), 
#              shape = 21, size = 1.5, stroke = 0.2, alpha = 0.6) +
#   geom_hline(yintercept=100, linetype='dashed', col = 'lightgrey')+
#   labs(title = "",
#        x = "Dataset", y = "% fragments mapping to correct gene orientation") +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
#         axis.title = element_text(size = 12, face = "bold"),
#         plot.title = element_text(size = 14, face = "bold"),
#         legend.position = "none") +
#   scale_x_discrete(labels = datasetsLabels) +
#   scale_fill_manual(values = datasetsPalette, labels = datasetsLabels) +
#   scale_color_manual(values = datasetsPalette, labels = datasetsLabels)
# 
# setwd("~/cfRNA-meta/slide_deck/")
# #setwd("~/ESMO")
# ggsave("mapping_expected_strand_with_flomics.png", p, width = 8, height = 6)
# ggsave("mapping_expected_strand_with_flomics.pdf", p, width = 8, height = 6)
# 
# 
# for (i in 1:length(ggplot_objects)) {
#   col_name <- column_names[i]
#   output_file <- paste0(gsub(" ", "_", tolower(col_name)), "_test_plot_with_wang.png")
#   ggsave(output_file, ggplot_objects[[i]], width = 8, height = 6)
#   output_file <- paste0(gsub(" ", "_", tolower(col_name)), "_test_plot_with_wang.pdf")
#   ggsave(output_file, ggplot_objects[[i]], width = 8, height = 6)
# }
# 
# x_bw <- bw.nrd0(joint_table_final$exonic_reads_minus_spike_ins)
# y_bw <- bw.nrd0(joint_table_final$Genes_contributing_to_80._of_reads)
# 
# p <- ggplot(data = joint_table_final, aes(x = exonic_reads_minus_spike_ins, y = Genes_contributing_to_80._of_reads, color = dataset)) +
#   geom_point(size = 2, alpha = 0.7) + 
#   # geom_density_2d(
#   #   aes(color = dataset, group = dataset),
#   #   alpha = 0.3, 
#   #   size = 0.8,
#   #   contour_var = "ndensity", # or "count", "density"
#   #   h = c(20, 3)         # optional manual bandwidth
#   # )+
#   geom_density_2d(
#     data = subset(joint_table_final, dataset == "liquidx"),
#     aes(color = dataset),
#     alpha = 0.5, size = 0.8
#   ) +
#   labs(title = "",
#        x = "% of reads mapping to exons",
#        y = "# of Genes Contributing to 80% of Reads",
#        color = "")+
#   theme_minimal(base_size = 14) + 
#   theme(
#     plot.title = element_text(size = 16, face = "bold", hjust = 0.5), 
#     axis.title = element_text(size = 14, face = "bold"), 
#     axis.text = element_text(size = 12), 
#     legend.title = element_text(size = 12, face = "bold"),
#     legend.text = element_text(size = 10), 
#     legend.position = "right", 
#     panel.grid.major = element_line(color = "gray80"), 
#     panel.grid.minor = element_blank(),
#     plot.background = element_rect(
#       fill = "white",
#       colour = "white"
#     )
#   ) +
#   scale_color_manual(values = datasetsPalette, labels= datasetsLabels) 
# 
# ps <- p  +  scale_y_continuous(trans=log2_trans()) 
# ps
# 
# ggsave("diversity_scatterplot_flomics_only_no_sevilla.png", plot = ps, width = 12, height = 8, dpi = 300)
# ggsave("diversity_scatterplot_flomics_only_no_sevilla.pdf", plot = ps, width = 12, height = 8, dpi = 300)
# ##### 
# #Plot scatter plots genes contributing to % of reads vs % exonic per dataset
# unique_batches <- unique(table_filtered$sequencing_batch)
# 
# for (batch in unique_batches) {
#   obj_dataset <- table_filtered[table_filtered$sequencing_batch == batch, ]
#   
#   p <- ggplot(data = obj_dataset, aes(x = exonic_percentage , y = log(genes_contributing_to_80._of_reads))) +
#     geom_point() +
#     labs(title = paste("Diversity scatterplot for", batch),
#          x = "exonic percentage",
#          y = "log(# of genes contributing to 80% of reads)") +
#     theme_bw()
#   
#   file_name <- paste0("plot_", batch, ".png")
#   
#   ggsave(filename = file_name, plot = p, width = 8, height = 6)
# }
