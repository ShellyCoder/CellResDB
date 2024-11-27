# Clear the workspace and set options
rm(list = ls())
options(stringsAsFactors = FALSE)

# Load required libraries
library(dplyr)
library(ggplot2)

# Load tissue information and sample datasets
tissue_info <- read.table("/Users/luca/myWorkSpace/web/cellresponse/data/tissue_Search.txt", sep = "\t", header = TRUE)
demo_samples <- read.table("/Users/luca/myWorkSpace/web/cellresponse/data/samples.txt", sep = "\t", header = TRUE)

# Remove unnecessary columns from tissue_info
tissue_info$diseaseTypeRecorded <- NULL
tissue_info$pubmedID <- NULL

# Remove duplicate samples caused by shared 'pre' and select relevant columns
demo_samples <- demo_samples[!duplicated(demo_samples$Sample.ID), 
                             c("Sample.ID_RAW", "Sample.ID", "Database.ID", "Cell.counts", "Response", "Drug.type")]

# Merge samples with tissue information based on 'Database.ID'
demo_sample_tissue <- dplyr::left_join(demo_samples, tissue_info, by = c("Database.ID" = "datasetID"))

# Standardize some disease names for consistent visualization
demo_sample_tissue[demo_sample_tissue$standardDisease == "B-ALL", "standardDisease"] <- "ALL"
demo_sample_tissue[demo_sample_tissue$standardDisease == "FLT3-ITD AML", "standardDisease"] <- "AML"
demo_sample_tissue[demo_sample_tissue$standardDisease == "IgG lambda MM", "standardDisease"] <- "MM"
demo_sample_tissue[demo_sample_tissue$standardDisease == "CM", "standardDisease"] <- "MEL"

# Filter out untreated samples and select relevant columns
demo_sample_tissue$sample_index <- paste0("sample", 1:nrow(demo_sample_tissue))
demo_sample_tissue_R_NR <- demo_sample_tissue[demo_sample_tissue$Response != "Untreated", ]
demo_sample_tissue_R_NR <- demo_sample_tissue_R_NR[, c("sample_index", "standardDisease", "Drug.type", "Response")]

# Convert Response to binary labels (1 for "R", 0 for "NR")
demo_sample_tissue_R_NR <- demo_sample_tissue_R_NR %>%
  mutate(Response = ifelse(Response == "R", 1, 0))

# Calculate total samples and response rates for each drug type and cancer type
summary_df <- demo_sample_tissue_R_NR %>%
  group_by(Drug.type, standardDisease) %>%
  summarise(
    total_samples = n(),  # Total samples
    response_rate = sum(Response) / n(),  # Proportion of samples with a response
    .groups = "drop"
  )

# Print the summary table for verification
print(summary_df)

# Create a bubble plot to visualize the relationship between cancer types and treatment combinations
ggplot(summary_df, aes(x = standardDisease, y = Drug.type)) +
  geom_point(
    aes(size = response_rate, color = log2(total_samples)), 
    alpha = 0.7  # Adjust transparency for better visibility
  ) +
  scale_size_continuous(range = c(1, 8), name = "Response Rate") +  # Scale bubble size by response rate
  scale_color_gradient(low = "lightblue", high = "darkblue", name = "log2(Total Samples)") +  # Color gradient for total samples
  theme_minimal() +  # Minimal theme for clean appearance
  labs(
    title = "Cancer Types and Treatment Combinations",
    x = "Cancer Type",
    y = "Treatment Combination"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability



