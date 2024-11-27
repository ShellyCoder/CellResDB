# Clear the workspace and set options
rm(list = ls())
options(stringsAsFactors = FALSE)

# Load required libraries
library(dplyr)
library(ggtree)
library(tidyverse)
library(data.tree)
library(ape)

# Load datasets
tissue_info <- read.table("/Users/luca/myWorkSpace/web/cellresponse/data/tissue_Search.txt", 
                          sep = "\t", header = TRUE)
demo_samples <- read.table("/Users/luca/myWorkSpace/web/cellresponse/data/samples.txt", 
                           sep = "\t", header = TRUE)

# Remove unnecessary columns from the tissue_info dataset
tissue_info$diseaseTypeRecorded <- NULL
tissue_info$pubmedID <- NULL

# Remove duplicate samples caused by shared "pre" IDs and select relevant columns
demo_samples <- demo_samples[!duplicated(demo_samples$Sample.ID), 
                             c("Sample.ID_RAW", "Sample.ID", "Database.ID", "Cell.counts", "Response")]
table(is.na(demo_samples))  # Check for missing values

# Merge samples with tissue information based on "Database.ID"
demo_sample_tissue <- dplyr::left_join(demo_samples, tissue_info, 
                                       by = c("Database.ID" = "datasetID"))
table(is.na(demo_sample_tissue))  # Check for missing values after merging

# Adjust cancer subtype names for consistent display
demo_sample_tissue[demo_sample_tissue$standardDisease == "B-ALL", "standardDisease"] <- "ALL"
demo_sample_tissue[demo_sample_tissue$standardDisease == "FLT3-ITD AML", "standardDisease"] <- "AML"
demo_sample_tissue[demo_sample_tissue$standardDisease == "IgG lambda MM", "standardDisease"] <- "MM"
demo_sample_tissue[demo_sample_tissue$standardDisease == "CM", "standardDisease"] <- "MEL"

# Inspect categorical distributions
table(demo_sample_tissue$Response)
table(demo_sample_tissue$standardDisease)
table(demo_sample_tissue$tissue)

# Create unique identifiers for each sample
demo_sample_tissue$sample_index <- paste0("sample", 1:nrow(demo_sample_tissue))

# Add placeholder columns for padding (used for visualization)
demo_sample_tissue$padding1 <- "padding1"
demo_sample_tissue$padding2 <- "padding2"
demo_sample_tissue$padding3 <- "padding3"

# Calculate the total number of cells across all samples
sum(demo_sample_tissue$Cell.counts)

# Arrange the data by tissue type, disease type, and response
demo_sample_tissue <- demo_sample_tissue %>%
  arrange(tissue, standardDisease, Response)

# Generate a star tree for visualizing the relationships
n_samples <- nrow(demo_sample_tissue)
tree <- stree(n_samples, type = "star")
tree$tip.label <- demo_sample_tissue$sample_index  # Assign sample indices as tip labels

# Extract attributes for samples to map colors in the plot
df <- demo_sample_tissue[, c("Response", "padding1", "standardDisease", 
                             "padding2", "padding3", "tissue")] %>%
  magrittr::set_rownames(demo_sample_tissue$sample_index)

# Generate a circular dendrogram using ggtree
circ <- ggtree(tree, layout = "circular")

# Define color mappings for different attributes
color_mapping <- c(
  "padding1" = "white", "padding2" = "white", "padding3" = "white",
  "Bladder" = "#ff7bab", "BLCA" = "#fcb6d0",
  "Blood" = "#ff7b7b", "ALL" = "#ffbbbb", "AML" = "#e29696", "CLL" = "#cc7575",
  "Bone" = "#fca67c", "MM" = "#e8b297",
  "Breast" = "#fcd37a", "TNBC" = "#f7dea6",
  "Colorectum" = "#cacc64", "CRC" = "#d1d38d", "DC" = "#f0f2af",
  "Kidney" = "#669966", "ccRCC" = "#64b264", "PRCC" = "#7fc67f",
  "Liver" = "#2a965e", "HCC" = "#4f9e74",
  "Lung" = "#68a097", "LUAD" = "#80b7ae", "LUSC" = "#75dbca",
  "Lymph node" = "#3a767d", "LBCL" = "#548b8e", "MALT" = "#5e8284",
  "Pelvic cavity" = "#336a7f", "EC" = "#3f839e",
  "Prostate" = "#495a91", "CSPC" = "#596eb2",
  "Skin" = "#18297c", "BCC" = "#283da5", "MEL" = "#3d52b7", "OM" = "#4c62cc", "SCC" = "#4a64e5",
  "NR" = "#dd3b33", "R" = "#006b9d", "Untreated" = "#c6b430"
)

# Add a heatmap to the dendrogram with sample attributes
gheatmap(circ, df, offset = 0.8, width = 1.5,
               colnames_angle = 95, colnames_offset_y = 0.25) +
  scale_fill_manual(
    values = color_mapping,
    name = "Sample Attributes"
  )
  





