# Clear the workspace and set options
rm(list = ls())
options(stringsAsFactors = FALSE)

# Load required libraries
library(dplyr)
library(ggplot2)

# Load the dataset
demo_dataset <- read.table("/Users/luca/myWorkSpace/web/cellresponse/data/demo_datasetInfo.txt", 
                           sep = "\t", header = TRUE)

# Select relevant columns for analysis
demo_dataset <- demo_dataset[, c("Database_ID", "Therapeutic_regimen", "Drug_type", "Cell.counts", "Sample.counts")]

# Check for missing values in the dataset
table(is.na(demo_dataset))

# Process the 'Therapeutic_regimen' column to count individual drug usage
drug_usage <- demo_dataset %>%
  # Split the therapeutic regimen string into individual drugs
  mutate(Therapeutic_regimen = strsplit(as.character(Therapeutic_regimen), " \\+ ")) %>%
  # Expand the list of drugs into individual rows
  unnest(Therapeutic_regimen) %>%
  # Count the frequency of each drug
  count(Therapeutic_regimen, sort = TRUE)

# Sort the drug usage data frame by usage count in descending order
drug_usage <- drug_usage[order(drug_usage$n, decreasing = TRUE), ] %>% as.data.frame()

# Convert the 'Therapeutic_regimen' column to a factor to maintain the order
drug_usage$Therapeutic_regimen <- factor(drug_usage$Therapeutic_regimen, levels = drug_usage$Therapeutic_regimen)

# Create a horizontal bar plot to visualize the distribution of drug usage
ggplot(drug_usage, aes(x = Therapeutic_regimen, y = n, fill = Therapeutic_regimen)) +
  geom_bar(stat = "identity") +  # Create bars based on drug usage counts
  geom_text(aes(label = n), vjust = -0.5, size = 3.5) +  # Add text labels above the bars
  coord_flip() +  # Flip coordinates to make a horizontal bar chart
  theme_minimal() +  # Use a minimal theme for the plot
  theme(
    legend.position = "none",  # Remove the legend
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels for better readability
  ) +
  labs(
    x = "Drug", 
    y = "Usage Count", 
    title = "Distribution of Drug Usage"
  )
