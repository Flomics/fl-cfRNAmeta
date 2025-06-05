setwd("~/cfRNA-meta/")

data <- read.table("../tables/cfRNA-meta_per_sample_metadata.tsv", header = TRUE, sep = "\t")
data_2nd_gen <- read.table("~/cfRNA-meta/sampleinfo_snakeDA_flomics_1.tsv", header = TRUE, sep = "\t")
data_2nd_gen$status <- "healthy"
data_liquidx <- read.table("~/cfRNA-meta/sampleinfo_snakeDA_liquidx_2025_05_08.tsv", header = TRUE, sep = "\t", fileEncoding = "UTF-8")

# Remove inserm CRC
data_liquidx <- data_liquidx[data_liquidx$collection_subcenter_health_care_unit != "INSERM-CRC23",]
#Remove "old" plasma sample age samples
data_liquidx <- data_liquidx[data_liquidx$plasma_sample_age < 4600,]
# Remove sevilla
data_liquidx <- data_liquidx[data_liquidx$collection_subcenter != "Sevilla",]
write.table(data_liquidx, "liquidx_filter_plasma_sample_age_no_sevilla.tsv", sep = "\t")
# Filter using selected 50 healthy samples
healthy_samples <- read.table("~/fl-cfRNAmeta/tables/healthy_matched.txt")
data_liquidx <- data_liquidx[data_liquidx$sample_name %in% healthy_samples[[1]], ]
#Remove non-healthy samples (for meta-analysis only)
#data_liquidx <- data_liquidx[data_liquidx$status == "healthy"]
# Remove samples not passing QC
#data_liquidx <- data_liquidx[data_liquidx$percentage_of_spliced_reads > 20,]

### Merge the two dataframes
merged_df_1 <- merge(data, data_2nd_gen, all = TRUE)
merged_df <- merge(merged_df_1, data_liquidx, all = TRUE)


### Merge the two dataframes
data_barplot <- merged_df
data_barplot$sequencing_batch[grepl("^FL", data_barplot$sequencing_batch)] <- "Flomics_2"

data_barplot$phenotype[grep("Acute Myeloid Leukemia",data_barplot$phenotype)] <- "cancer"
data_barplot$phenotype[grep("Chronic hepatitis B",data_barplot$phenotype)] <- "non-cancer disease"
data_barplot$phenotype[grep("Chronic kidney failure EPO-treated",data_barplot$phenotype)] <- "non-cancer disease"
data_barplot$phenotype[grep("Cirrhosis",data_barplot$phenotype)] <- "non-cancer disease"
data_barplot$phenotype[grep("colorectal cancer",data_barplot$phenotype)] <- "cancer"
data_barplot$phenotype[grep("Control",data_barplot$phenotype)] <- "healthy"
data_barplot$phenotype[grep("Diffuse large B-cell lymphoma",data_barplot$phenotype)] <- "cancer"
data_barplot$phenotype[grep("Esophagus cancer",data_barplot$phenotype)] <- "cancer"
data_barplot$phenotype[grep("G-CSF-treated healthy donors",data_barplot$phenotype)] <- "healthy"
data_barplot$phenotype[grep("Liver cancer",data_barplot$phenotype)] <- "cancer"
data_barplot$phenotype[grep("Lung cancer",data_barplot$phenotype)] <- "cancer"
data_barplot$phenotype[grep("Multiple myeloma",data_barplot$phenotype)] <- "cancer"
data_barplot$phenotype[grep("Nonalcoholic fatty liver disease",data_barplot$phenotype)] <- "non-cancer disease"
data_barplot$phenotype[grep("Nonalcoholic steatohepatitis",data_barplot$phenotype)] <- "non-cancer disease"
data_barplot$phenotype[grep("Pancreatic cancer",data_barplot$phenotype)] <- "cancer"
data_barplot$phenotype[grep("Pre-eclampsia",data_barplot$phenotype)] <- "non-cancer disease"
data_barplot$phenotype[grep("Primary mediastinal B-cell lymphoma",data_barplot$phenotype)] <- "cancer"
data_barplot$phenotype[grep("Stomach cancer",data_barplot$phenotype)] <- "cancer"
#data_barplot$status[is.na(data_barplot$status)] <- "missing"
#data_barplot$phenotype[is.na(data_barplot$phenotype)] <- "missing"

# Merge phenotype with status when phenotype is NA and status is not
data_barplot$phenotype[is.na(data_barplot$phenotype) & !is.na(data_barplot$status)] <- 
  data_barplot$status[is.na(data_barplot$phenotype) & !is.na(data_barplot$status)]

# Replace remaining NA values in phenotype with "missing"
data_barplot$phenotype[is.na(data_barplot$phenotype)] <- "missing"
# Replace empty strings in phenotype with "missing"
data_barplot$phenotype[data_barplot$phenotype == ""] <- "missing"



table(data_barplot$phenotype)
# Replace sequencing_batch with dataset_short_name when dataset_short_name is not NA
data_barplot$sequencing_batch[!is.na(data_barplot$dataset_short_name)] <-
  data_barplot$dataset_short_name[!is.na(data_barplot$dataset_short_name)]


plot_data <- data_barplot %>%
  group_by(sequencing_batch, phenotype) %>%
  summarise(count = n(), .groups = 'drop')

library(RColorBrewer)

custom_colors <- brewer.pal(4, "Dark2")
print(custom_colors)


ggplot(plot_data, aes(x = sequencing_batch, y = count, fill = phenotype)) +
  geom_bar(stat = "identity") +
  labs(title = "Sample Status per Dataset", x = "Dataset", y = "Count") +
  scale_fill_manual(values = brewer.pal(4, "Dark2")) +
  theme_minimal(base_size = 20) +  # Increase base text size for better readability
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 18, color = "black"),
    axis.title.x = element_text(size = 20, face = "bold"), 
    axis.title.y = element_text(size = 20, face = "bold"),  
    legend.title = element_text(size = 18, face = "bold"),  
    legend.text = element_text(size = 16, face = "bold"),  
    panel.grid.major = element_line(size = 0.8),  # Thicker major grid lines for clarity
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    plot.margin = margin(20, 20, 20, 20)  # Increase plot margins for better spacing
  ) +
  labs(title = "",
       x = "Dataset",
       y = "Number of samples")

