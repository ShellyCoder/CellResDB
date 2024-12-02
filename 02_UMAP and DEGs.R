# Clear the environment
rm(list = ls())
options(stringsAsFactors = F)

# Load necessary libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(clusterProfiler)  
library(org.Hs.eg.db)     
library(dplyr)     

setwd("/sc/processed_data/")
file_list <- list.files(pattern = "rds$")

# Create output directory for UMAP plots
umap_output_dir <- "umap"
dir.create(umap_output_dir)

# Function to calculate legend text size dynamically based on the number of clusters
calculate_legend_text_size <- function(num_clusters) {
  base_size <- 12  # Base text size
  base_clusters <- 7  # Base number of clusters
  
  # Inversely scale the text size based on the number of clusters
  scale_factor <- base_clusters / num_clusters
  text_size <- base_size * scale_factor
  
  min_size <- 6  # Minimum text size
  max_size <- 12  # Maximum text size
  text_size <- max(min_size, min(text_size, max_size))
  
  return(text_size)
}

# Function to extract the legend from a ggplot object
get_legend <- function(myggplot) {
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Loop over each file in the file_list
for (file in file_list) {
  # Load Seurat object
  data <- readRDS(file)
  
  # Set cluster identity for the Seurat object
  Idents(data) <- data@meta.data$cluster
  
  # Calculate the number of clusters
  num_clusters <- length(unique(data$cluster))
  
  # Dynamically calculate legend text size based on the number of clusters
  legend_text_size <- calculate_legend_text_size(num_clusters)
  
  # Run UMAP on the Seurat object
  data <- RunUMAP(data, dims = 1:20, verbose = TRUE)
  
  # Set a common theme for the plots
  common_theme <- theme(
    plot.title = element_blank(),  # Remove plot title
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")  # Set uniform margins
  )
  
  # UMAP plot colored by cluster
  p1 <- DimPlot(data, reduction = "umap", group.by = "cluster") +
    common_theme +
    theme(legend.position = "right", legend.text = element_text(size = 20)) +  # Fixed legend position and size
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))  # Vertical legend
  
  # Extract the legend from the plot
  p1_legend <- get_legend(p1)
  p1 <- p1 + theme(legend.position = "none")  # Remove legend from the plot
  
  # Combine UMAP plot and legend, then save as .jpg
  umap_plot <- plot_grid(p1, p1_legend, rel_widths = c(2, 1))
  print(umap_plot)
  ggsave(paste0(umap_output_dir, "/umap_cluster_", file_path_sans_ext(file), ".jpg"), plot = umap_plot, width = 20, height = 10)
  
  # UMAP plot colored by treatment
  p2 <- DimPlot(data, reduction = "umap", group.by = "treatment") +
    common_theme +
    theme(legend.position = "right", legend.text = element_text(size = 18))  # Treatment legend
  
  # Extract the legend from the second plot
  p2_legend <- get_legend(p2)
  p2 <- p2 + theme(legend.position = "none")  # Remove legend from the plot
  
  # Combine UMAP plot and legend, then save as .jpg
  umap_plot <- plot_grid(p2, p2_legend, rel_widths = c(2, 1))
  print(umap_plot)
  ggsave(paste0(umap_output_dir, "/umap_treatment_", file_path_sans_ext(file), ".jpg"), plot = umap_plot, width = 10, height = 5)
  
  # Clean up memory
  rm(data)
  gc()
}

# DEGs
# Directory for saving the marker results
allmarker_output_dir <- "allmarker"
dir.create(allmarker_output_dir)

# Loop over each file in the file list
for (file in file_list) {
  tryCatch({
    # Load Seurat object
    data <- readRDS(file)
    
    # Set the cluster identity
    Idents(data) <- "cluster"
    
    # Find all markers for each cluster
    markers <- FindAllMarkers(data, test.use = "wilcox", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    
    # Filter markers with adjusted p-value < 0.05 and select gene column first
    all.markers <- markers %>%
      dplyr::select(gene, everything()) %>%
      subset(p_val_adj < 0.05)
    
    # Convert gene symbols to ENTREZ IDs
    ids <- bitr(all.markers$gene, 'SYMBOL', 'ENTREZID', 'org.Hs.eg.db')
    
    # Merge the ENTREZ ID data with the marker data, genes without ENTREZ ID will be filtered out
    all.markers <- merge(all.markers, ids, by.x = 'gene', by.y = 'SYMBOL')
    
    # Sort markers by cluster
    sorted_markers <- all.markers %>%
      arrange(cluster)
    
    # Save the sorted markers to a CSV file
    write.csv(sorted_markers, paste0(allmarker_output_dir, "/allmarkers_", file_path_sans_ext(file), ".csv"), row.names = TRUE)
    
    # Clean up memory
    rm(data)
    gc()
  }, error = function(e) {
    # Catch any errors and print the file name and error message
    print(paste0(file, ": ", e))
  })
}






