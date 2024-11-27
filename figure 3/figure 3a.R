# Clear the workspace
rm(list = ls())
options(stringsAsFactors = FALSE)

# Load required libraries
library(dplyr)
library(ggtree)
library(tidyverse)
library(data.tree)
library(ape)

# Load datasets
demo_dataset <- read.table("/Users/luca/myWorkSpace/web/cellresponse/data/demo_datasetInfo.txt", 
                           sep = "\t", header = TRUE)
tissue_info <- read.table("/Users/luca/myWorkSpace/web/cellresponse/data/tissue_Search.txt", 
                          sep = "\t", header = TRUE)

# Remove unnecessary columns from the tissue_info dataset
tissue_info$diseaseTypeRecorded <- NULL
tissue_info$pubmedID <- NULL

# Select specific columns from the demo_dataset for analysis
demo_dataset <- demo_dataset[, c("Database_ID", "Therapeutic_regimen", 
                                 "Drug_type", "Cell.counts", "Sample.counts")]
table(is.na(demo_dataset))  # Check for missing values

# Merge datasets based on the common column 'Database_ID'
demo_dataset_tissue <- dplyr::left_join(demo_dataset, tissue_info, 
                                        by = c("Database_ID" = "datasetID"))
table(is.na(demo_dataset_tissue))  # Check for missing values after merging

# Adjust cancer subtype names for better visualization
demo_dataset_tissue[demo_dataset_tissue$standardDisease == "B-ALL", "standardDisease"] <- "ALL"
demo_dataset_tissue[demo_dataset_tissue$standardDisease == "FLT3-ITD AML", "standardDisease"] <- "AML"
demo_dataset_tissue[demo_dataset_tissue$standardDisease == "IgG lambda MM", "standardDisease"] <- "MM"
demo_dataset_tissue[demo_dataset_tissue$standardDisease == "CM", "standardDisease"] <- "MEL"

# Add placeholder columns for padding (used for visualization purposes)
demo_dataset_tissue$padding1 <- "padding1"
demo_dataset_tissue$padding2 <- "padding2"

# Prepare data for hierarchical clustering
data_num <- data.frame(
  standardDisease = as.numeric(factor(demo_dataset_tissue$standardDisease)),
  tissue = as.numeric(factor(demo_dataset_tissue$tissue)),
  padding1 = as.numeric(factor(demo_dataset_tissue$padding1))
)

rownames(data_num) <- demo_dataset_tissue$Database_ID

# Compute distance matrix using Euclidean distance
dist_matrix <- dist(data_num, method = "euclidean")

# Perform hierarchical clustering
hclust_res <- hclust(dist_matrix, method = "ward.D2")

# Convert clustering results to a phylogenetic tree
tree <- as.phylo(hclust_res)

# Extract attributes for samples to map colors in the plot
df <- demo_dataset_tissue[, c("standardDisease", "padding1", "tissue")] %>%
  magrittr::set_rownames(demo_dataset_tissue$Database_ID)

# Generate a circular dendrogram using ggtree
circ <- ggtree(tree, layout = "circular")

# Define color mappings for different attributes
color_mapping <- c(
  "padding1" = "white", "padding2" = "white",
  "Bladder" = "#ff7bab", "BLCA" = "#fcb6d0",
  "Blood" = "#ff7b7b", "ALL" = "#ffbbbb", "AML" = "#e29696", "CLL" = "#cc7575",
  "Bone" = "#fca67c", "MM" = "#e8b297",
  "Breast" = "#fcd37a",  "TNBC" = "#f7dea6",
  "Colorectum" = "#cacc64",  "CRC" = "#d1d38d",  "DC" = "#f0f2af",
  "Kidney" = "#669966",  "ccRCC" = "#64b264",  "PRCC" = "#7fc67f",
  "Liver" = "#2a965e",  "HCC" = "#4f9e74",
  "Lung" = "#68a097",  "LUAD" = "#80b7ae",  "LUSC" = "#75dbca",
  "Lymph node" = "#3a767d",  "LBCL" = "#548b8e",  "MALT" = "#5e8284",
  "Pelvic cavity" = "#336a7f", "EC" = "#3f839e",
  "Prostate" = "#495a91", "CSPC" = "#596eb2",
  "Skin" = "#18297c", "BCC" = "#283da5", "MEL" = "#3d52b7", "OM" = "#4c62cc", "SCC" = "#4a64e5",
  "yes" = "#a1a5b5"
)

# Add heatmap to the dendrogram with sample attributes
gheatmap(circ, df, offset = 0.8, width = 1.8, 
         colnames_angle = 95, colnames_offset_y = 0.25) +
  scale_fill_manual(
    values = color_mapping,
    name = "Sample Attributes"
  )
