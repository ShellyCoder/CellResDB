# Clear the environment
rm(list = ls())
options(stringsAsFactors = F)

library(Seurat)
library(RColorBrewer)
library(dplyr)
library(magrittr)
library(CellChat)
library(patchwork)
library(tidydr)
library(igraph)
library(ggrepel)
library(pheatmap)

setwd("sc/processed_data/")
file_list <- list.files(pattern = "rds$")
dir.create("cellchat")

for (file in file_list) {
  tryCatch({
    # Read the RDS file
    data <- readRDS(file)
    
    # Get all unique treatment groups
    groups <- unique(data@meta.data$treatment)
    
    # Extract file name without extension
    file_name <- gsub("(.*)\\.rds$", "\\1", file)
    output_dir <- paste0("cellchat/", file_name)
    
    # Create output directory if it does not exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir)
    }
    
    # Loop over each treatment group
    for (group in groups) {
      # Get cell IDs for the current treatment group
      cell_ids <- rownames(data@meta.data[data@meta.data$treatment == group, ])
      
      # Subset the data for the current treatment group
      R_object <- subset(data, cells = cell_ids)
      
      # Extract the RNA assay data
      R.data.input <- GetAssayData(R_object, assay = "RNA", slot = "data")
      
      # Extract relevant metadata (cluster and treatment)
      R.meta <- R_object@meta.data[, c("cluster", "treatment")]
      
      # Create the CellChat object using the expression data
      R.cellchat <- createCellChat(object = R.data.input)
      
      # Add metadata to the CellChat object
      R.cellchat <- addMeta(R.cellchat, meta = R.meta)
      
      # Set the cluster information as the identity class
      R.cellchat <- setIdent(R.cellchat, ident.use = "cluster")
      
      # Set the database for human data (use mouse database if applicable)
      R.cellchat@DB <- CellChatDB.human  # For mouse, use CellChatDB.mouse
      
      # Subset the data (all genes are included by default)
      R.cellchat <- subsetData(R.cellchat, features = NULL)
      
      # Set the plan for parallel processing with one worker
      future::plan("multisession", workers = 1)
      
      # Increase the maximum memory size limit (optional, adjust as needed)
      options(future.globals.maxSize = 80 * 1024^3)  # Set to 80GB, adjust as needed
      
      # Identify over-expressed genes (ligands and receptors)
      R.cellchat <- identifyOverExpressedGenes(R.cellchat)
      
      # Identify interactions between over-expressed genes (pathways)
      R.cellchat <- identifyOverExpressedInteractions(R.cellchat)
      
      # Project gene data onto the Protein-Protein Interaction (PPI) network
      R.cellchat <- projectData(R.cellchat, PPI.human)  # For mouse, use PPI.mouse
      
      # Compute communication probabilities for the raw data
      R.cellchat <- computeCommunProb(R.cellchat, raw.use = TRUE)
      
      # Filter out networks with fewer than 10 cells (irrelevant interactions)
      R.cellchat <- filterCommunication(R.cellchat, min.cells = 10)
      
      # Compute communication probabilities at the pathway level
      R.cellchat <- computeCommunProbPathway(R.cellchat)
      
      # Aggregate the network
      R.cellchat <- aggregateNet(R.cellchat)
      
      # Compute centrality of the network
      R.cellchat <- netAnalysis_computeCentrality(R.cellchat, slot.name = "netP")
      
      # Save the processed CellChat object for the current group
      saveRDS(R.cellchat, file = paste0(output_dir, "/R_cellchat_", group, ".rds"))
    }
  }, error = function(e) {
    # Catch and print any errors
    cat("Error in file", file, ": ", e$message, "\n")
  })
}

setwd("cellchat/")
dir_list <- list.dirs()
dir_list <- dir_list[-1]

for (dir in dir_list) {
  setwd(dir)  # Set working directory to current dir in dir_list
  file_list <- list.files(pattern = "rds$")  # Get list of all .rds files
  
  # Create necessary directories if they don't exist
  if (!dir.exists("heatmap")) {
    dir.create("heatmap")  
  }
  
  for (file in file_list) {
    # Read the RDS file (CellChat object)
    R.cellchat <- readRDS(file)
    groupSize <- as.numeric(table(R.cellchat@idents))  # Cell group sizes
    group <- gsub("R_cellchat_(.*)\\.rds", "\\1", file)  # Extract group name from the filename
    
    # ----------- Part 1: Interaction Network Plot (JPEG) -----------
    # Create a JPEG device for the first plot (interaction network)
    jpeg(paste0("interaction_network_", group, "_count.jpg"), width = 8, height = 8, units = 'in', res = 300)
    
    # Plot the interaction network
    par(mar = c(5, 4, 4, 6) + 0.1, xpd = TRUE)
    netVisual_circle(R.cellchat@net$count, 
                     vertex.weight = groupSize, 
                     weight.scale = TRUE, 
                     label.edge = FALSE,
                     title.name = "Number of interactions",
                     vertex.label.cex = 1)
    
    # Close the JPEG device
    dev.off()
    
    # ----------- Part 2: Heatmap Plot (JPEG) -----------
    output_dir <- paste0("heatmap/", dir)
    
    # Create output directory for heatmaps if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir)
    }
    
    # Generate file name for the heatmap
    file_name <- paste0(output_dir, "/CellChat_Network_", group, "_heatmap.jpg")
    
    # Get the interaction count data (matrix)
    heatmap_data <- R.cellchat@net$count
    num_rows <- nrow(heatmap_data)
    num_cols <- ncol(heatmap_data)
    
    # Set dynamic cell size based on heatmap dimensions
    base_size <- 20
    cellwidth <- max(5, min(base_size, 1000 / num_cols))  # Column cell width
    cellheight <- max(5, min(base_size, 800 / num_rows))  # Row cell height
    
    # Create a JPEG device for the heatmap
    jpeg(file_name, width = 4, height = 4, res = 300, units = "in")
    
    # Plot the heatmap
    p <- pheatmap(heatmap_data,
                  clustering_distance_rows = "euclidean",
                  clustering_distance_cols = "euclidean",
                  cluster_rows = TRUE,
                  cluster_cols = TRUE,
                  show_rownames = TRUE,
                  show_colnames = TRUE,
                  fontsize = 10, 
                  fontsize_row = 8,
                  fontsize_col = 8,
                  cellwidth = cellwidth,
                  cellheight = cellheight,
                  display_numbers = TRUE, 
                  number_format = "%.2f", 
                  main = "Interaction Counts",
                  treeheight_row = 2.5, 
                  treeheight_col = 2.5)
    
    # Save and close the heatmap plot
    print(p)
    dev.off()
    
    # Print a message for each file processed
    cat("Processed and saved:", file, "\n")
  }
  
  # Return to the original directory after processing the current directory
  setwd("sc/processed_data/cellchat")
  
  # Print completion message for the directory
  print(paste0(dir, " processed and saved"))
}

