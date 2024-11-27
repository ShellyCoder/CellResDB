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

# Select relevant columns and calculate log10 transformation of cell counts
demo_dataset <- demo_dataset[, c("Database_ID", "Therapeutic_regimen", "Drug_type", "Cell.counts", "Sample.counts")]
demo_dataset$Cell.log10 <- log10(demo_dataset$Cell.counts)

# Check for missing values in the dataset
table(is.na(demo_dataset))

# Create bins based on log10-transformed cell counts
demo_dataset <- demo_dataset %>%
  mutate(Bin = cut(
    Cell.log10, 
    breaks = seq(0, max(Cell.log10) + (1 - max(Cell.log10) %% 1), by = 1),  # Create bins of size 1
    right = FALSE, 
    include.lowest = TRUE
  ))

# Summarize the number of files in each bin
bin_summary <- demo_dataset %>%
  group_by(Bin) %>%
  summarise(Cell_Count = n()) %>%
  ungroup()

# Create a bar plot to show the distribution of files by cell count bins
ggplot(bin_summary, aes(x = Bin, y = Cell_Count, fill = Bin)) +
  geom_bar(stat = "identity") +  # Create bars based on the file count
  geom_text(aes(label = Cell_Count), vjust = -0.5, size = 3.5) +  # Add text labels above bars
  theme_minimal() +  # Use a minimal theme for the plot
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for readability
    legend.position = "none"  # Remove the legend
  ) +
  labs(
    x = "Log10 Cell Count (Bins)", 
    y = "Number of Files", 
    title = "Distribution of Files by Cell Count"
  ) +
  ggsci::scale_fill_npg()  # Use Nature Publishing Group color palette
