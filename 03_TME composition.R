# Clear the environment
rm(list = ls())
options(stringsAsFactors = F)

setwd("/sc/processed_data/")
file_list <- list.files(pattern = "rds$")

library(ggplot2)
library(dplyr)
library(Seurat)
library(cellPCT)
library(tidyr)
library(cowplot)
library(tools)
library(cols4all)
library(ggalluvial)


# Create directories for output
ratio_output_dir <- "ratio"
dir.create(ratio_output_dir)
table <- "ratio/table"
dir.create(table)
image <- "ratio/image"
dir.create(image)

# Loop through each file in the file list
for (file in file_list) {
  tryCatch({
    # Read the Seurat object
    data <- readRDS(file)
    
    # Get the number of unique clusters
    type_n <- length(unique(data@meta.data$cluster))
    
    # Generate dynamic colors for the clusters
    cols <- c4a('dynamic', type_n)
    
    # Function to extract the legend from a ggplot object
    get_legend <- function(myggplot) {
      tmp <- ggplot_gtable(ggplot_build(myggplot))
      leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
      legend <- tmp$grobs[[leg]]
      return(legend)
    }
    
    # Extract metadata from the Seurat object
    metadata <- data@meta.data
    
    # Calculate the number of cells in each sample-cluster combination
    cell_counts <- table(metadata$sample_ID, metadata$cluster)
    
    # Calculate the proportion of cells in each cluster for each sample
    cell_ratios <- prop.table(cell_counts, 1)  # Calculate proportions by row (sample)
    
    # Convert the proportions to a data frame
    cell_ratios_df <- as.data.frame(cell_ratios)
    
    # Add treatment and response information
    cell_ratios_df$treatment <- metadata$treatment[match(cell_ratios_df$Var1, metadata$sample_ID)]
    cell_ratios_df$response <- metadata$Response[match(cell_ratios_df$Var1, metadata$sample_ID)]
    cell_ratios_df$treatment_response <- paste(cell_ratios_df$treatment, cell_ratios_df$response, sep = "-")
    
    # Factor treatment and response columns for proper ordering
    cell_ratios_df$treatment <- factor(cell_ratios_df$treatment, levels = c("Pre", "Post"))
    cell_ratios_df$treatment_response <- factor(cell_ratios_df$treatment_response, levels = c("Pre-R", "Pre-NR", "Post-R", "Post-NR"))
    
    # Rename columns for clarity
    colnames(cell_ratios_df) <- c("sample", "cluster", "ratio", "treatment", "response", "treatment_response")
    
    # Sort by treatment response
    cell_ratios_df <- cell_ratios_df %>%
      arrange(treatment_response)
    
    # Save the cell ratios as a CSV file
    write.csv(cell_ratios_df, paste0(table, "/cell_ratio_", file_path_sans_ext(file), ".csv"))
    
    # Pivot the data to wide format
    wide_ratio_sorted <- pivot_wider(data = cell_ratios_df, names_from = sample, values_from = ratio, id_cols = cluster)
    wide_ratio_sorted <- wide_ratio_sorted %>%
      dplyr::rename(plot_cell_type = cluster)
    
    wide_ratio_sorted <- as.data.frame(wide_ratio_sorted)
    
    # Create Sankey plot based on the number of columns (samples)
    if (length(unique(colnames(wide_ratio_sorted))) - 1 < 8) {
      p1 <- plot_sankey(wide_ratio_sorted, coord_flip = TRUE) + 
        theme(
          text = element_text(size = 22),  # Global text size
          axis.title = element_text(size = 30),  # Axis title size
          axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust = 1),  # X-axis text size and angle
          axis.text.y = element_text(size = 20),  # Y-axis text size
          legend.text = element_text(size = 40)   # Legend text size
        ) + 
        guides(fill = guide_legend(ncol = 1))
      
      p1_legend <- get_legend(p1)
      p1 <- p1 + theme(legend.position = "none")
      
      # Combine the plot and legend, then save as a JPEG
      ratio_plot <- plot_grid(p1, p1_legend, rel_widths = c(3, 1))
      ggsave(paste0(image, "/", file_path_sans_ext(file), ".jpg"), plot = ratio_plot, width = 35, height = 15, dpi = 300, device = "jpeg", limitsize = FALSE)
    } else {
      p1 <- plot_sankey(wide_ratio_sorted) + 
        theme(
          text = element_text(size = 30),  # Global text size
          axis.title = element_text(size = 30),  # Axis title size
          axis.text.x = element_text(size = 5, angle = 45, vjust = 1, hjust = 1),  # X-axis text size and angle
          axis.text.y = element_text(size = 30),  # Y-axis text size
          legend.text = element_text(size = 66)  # Legend text size
        ) + 
        guides(fill = guide_legend(ncol = 1))
      
      p1_legend <- get_legend(p1)
      p1 <- p1 + theme(legend.position = "none")
      
      # Combine the plot and legend, then save as a JPEG
      ratio_plot <- plot_grid(p1, p1_legend, rel_widths = c(3, 1))
      ggsave(paste0(image, "/", file_path_sans_ext(file), ".jpg"), plot = ratio_plot, width = 35, height = 15, dpi = 300, device = "jpeg", limitsize = FALSE)
    }
    
    # Clean up memory after processing each file
    rm(data)
    gc()
    
  }, error = function(e) {
    # Handle errors by printing the file name and error message
    print(paste0(file, ": ", e))
  })
}
