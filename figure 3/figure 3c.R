# Clear the workspace and set options
rm(list = ls())
options(stringsAsFactors = FALSE)

# Load required libraries
library(dplyr)
library(ggplot2)
library(ggsci)  # For additional color palettes

# Load the dataset
demo_dataset <- read.table("/Users/luca/myWorkSpace/web/cellresponse/data/demo_datasetInfo.txt", 
                           sep = "\t", header = TRUE)

# Select relevant columns for analysis
demo_dataset <- demo_dataset[, c("Database_ID", "Therapeutic_regimen", "Drug_type", "Cell.counts", "Sample.counts")]

# Check for missing values in the dataset
table(is.na(demo_dataset))

# Process the 'Drug_type' column
# Split drug combinations, sort the components alphabetically, and recombine them
demo_dataset <- demo_dataset %>%
  mutate(Drug_type_sorted = sapply(strsplit(as.character(Drug_type), " \\+ "), function(x) {
    paste(sort(x), collapse = " + ")  # Sort and recombine drug names
  }))

# Count the frequency of each unique drug combination
drug_type_table <- table(demo_dataset$Drug_type_sorted)
drug_type_df <- as.data.frame(drug_type_table)
colnames(drug_type_df) <- c("Drug_type", "Count")  # Rename columns for clarity

# Sort the data frame by drug type for better visualization
drug_type_df <- drug_type_df[order(drug_type_df$Drug_type), ]

# Plotting the bar chart
ggplot(drug_type_df, aes(x = Drug_type, y = Count, fill = Drug_type)) +
  geom_bar(stat = "identity") +  # Create bar chart
  geom_text(aes(label = Count), vjust = -0.5, size = 3.5) +  # Add text labels above bars
  # coord_flip() +  # Uncomment to use horizontal bars for long labels
  theme_minimal() +  # Use a minimal theme for the plot
  theme(
    legend.position = "none",  # Remove legend
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels for readability
  ) +
  labs(
    x = "Drug Type", 
    y = "Count", 
    title = "Frequency of Drug Types"
  ) +
  ggsci::scale_fill_futurama()  # Use a Futurama color palette for bar fills

