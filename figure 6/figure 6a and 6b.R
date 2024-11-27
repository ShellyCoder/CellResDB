# Clear the workspace and set options
rm(list = ls())
options(stringsAsFactors = FALSE)

# Load required libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggtree)
library(aplot)
library(data.tree)
library(ape)

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

# Standardize some disease names for better visualization
demo_sample_tissue[demo_sample_tissue$standardDisease == "B-ALL", "standardDisease"] <- "ALL"
demo_sample_tissue[demo_sample_tissue$standardDisease == "FLT3-ITD AML", "standardDisease"] <- "AML"
demo_sample_tissue[demo_sample_tissue$standardDisease == "IgG lambda MM", "standardDisease"] <- "MM"
demo_sample_tissue[demo_sample_tissue$standardDisease == "CM", "standardDisease"] <- "MEL"

# Create a subset excluding untreated samples
demo_sample_tissue$sample_index <- paste0("sample", 1:nrow(demo_sample_tissue))
demo_sample_tissue_R_NR <- demo_sample_tissue[demo_sample_tissue$Response != "Untreated", ]

# Calculate counts and relative proportions of responses by disease type
response_counts <- demo_sample_tissue_R_NR %>%
  group_by(standardDisease, Response) %>%
  summarise(count = n(), .groups = "drop")

relative_proportions <- response_counts %>%
  group_by(standardDisease) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

# Plot relative proportions of NR and R responses (Figure 3A)
ggplot(relative_proportions, aes(x = standardDisease, y = proportion, fill = Response)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "Relative Proportion of NR and R Responses by Disease Type",
       x = "Disease Type", y = "Proportion", fill = "Response") +
  theme_minimal()

# Hierarchical clustering based on NR and R proportions
relative_df <- response_counts %>%
  pivot_wider(names_from = Response, values_from = count, values_fill = 0) %>%
  mutate(
    NR_proportion = NR / (NR + R),
    R_proportion = R / (NR + R)
  ) %>%
  as.data.frame()

# Create a distance matrix and cluster tree
data_num <- relative_df[, c("NR_proportion", "R_proportion")]
rownames(data_num) <- relative_df$standardDisease
dist_matrix <- dist(data_num, method = "euclidean")
hclust_res <- hclust(dist_matrix, method = "ward.D2")
tree <- as.phylo(hclust_res)

# Prepare data for heatmap
df <- demo_sample_tissue_R_NR[, c("standardDisease", "tissue")] %>%
  unique() %>%
  magrittr::set_rownames(demo_sample_tissue_R_NR$standardDisease)

# Create a rectangular hierarchical tree plot
circ <- ggtree(tree, layout = "rectangular")

# Define color mapping
color_mapping <- c(
  "Bladder" = "#ff7bab", "BLCA" = "#fcb6d0", "Blood" = "#ff7b7b", 
  "ALL" = "#ffbbbb", "AML" = "#e29696", "CLL" = "#cc7575", "Bone" = "#fca67c", 
  "MM" = "#e8b297", "Breast" = "#fcd37a", "TNBC" = "#f7dea6", "Colorectum" = "#cacc64", 
  "CRC" = "#d1d38d", "Kidney" = "#669966", "ccRCC" = "#64b264", 
  "Prostate" = "#495a91", "Skin" = "#18297c", "MEL" = "#3d52b7"
)

# Add heatmap to the tree
gheatmap(circ, df, offset = 0.8, width = 1.2,
               colnames_angle = 95, colnames_offset_y = 0.25) +
  scale_fill_manual(values = color_mapping, name = "Sample Attributes")


