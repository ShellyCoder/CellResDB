# Clear the workspace and set options
rm(list = ls())
options(stringsAsFactors = FALSE)

# Load required libraries
library(dplyr)
library(ggplot2)
library(ggsci)  # For additional color palettes

# Define the directory containing marker gene files
targeted_dir <- "/Users/luca/myWorkSpace/web/cellresponse/data/from_qiao/marker_gene/"
file_list <- list.files(path = targeted_dir, pattern = "\\.csv$", full.names = TRUE)

# Initialize a data frame to store the summary results
summary_df <- data.frame(File = character(), Cell_Type_Count = integer(), stringsAsFactors = FALSE)

# Loop through each CSV file
for (file in file_list) {
  # Read the data from the file
  data <- read.csv(file)
  
  # Count the number of unique cell types (based on 'cluster')
  cell_type_count <- data %>%
    distinct(cluster) %>%
    nrow()
  
  # Append the result to the summary data frame
  summary_df <- rbind(summary_df, data.frame(File = basename(file), Cell_Type_Count = cell_type_count))
}

# Bin the cell type counts into intervals of size 5 for grouping
summary_df <- summary_df %>%
  mutate(Bin = cut(
    Cell_Type_Count, 
    breaks = seq(0, max(Cell_Type_Count) + (5 - max(Cell_Type_Count) %% 5), by = 5), 
    right = FALSE, include.lowest = TRUE
  ))

# Summarize the number of files in each bin
bin_summary <- summary_df %>%
  group_by(Bin) %>%
  summarise(File_Count = n()) %>%
  ungroup()

# Create a bar plot to show the distribution of files by cell type count
ggplot(bin_summary, aes(x = Bin, y = File_Count, fill = Bin)) +
  geom_bar(stat = "identity") +  # Create bars based on the file count
  geom_text(aes(label = File_Count), vjust = -0.5, size = 3.5) +  # Add text labels above bars
  theme_minimal() +  # Use a minimal theme for the plot
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for readability
    legend.position = "none"  # Remove the legend
  ) +
  labs(
    x = "Number of Cell Types (Bins)", 
    y = "Number of Files", 
    title = "Distribution of Files by Cell Type Count"
  ) +
  ggsci::scale_fill_npg()  # Use Nature Publishing Group color palette
