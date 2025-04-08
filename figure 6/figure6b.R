# This script processes and visualizes the enrichment results (odds ratio and adjusted p-values)
# of different cell types across multiple datasets. It summarizes enrichment patterns and 
# plots a significance vs. dataset coverage bubble chart.

# Prepare environment
clean_result <- final_results[!sapply(final_results, is.null)]
related_df <- do.call("rbind", clean_result)
related_df <- dplyr::left_join(related_df, datasetInfo, by = "datasetID")

library(tidyr)
library(dplyr)

# Reshape adjusted p-values to wide format (cell types as rows, datasets as columns)
pval_selected <- related_df %>%
  dplyr::select(datasetID, cluster1, p_adj)

pval_wide <- pval_selected %>%
  pivot_wider(names_from = datasetID, values_from = p_adj) %>%
  as.data.frame()

rownames(pval_wide) <- pval_wide$cluster1
pval_wide$cluster1 <- NULL

# Reshape odds ratio values to wide format
OR_selected <- related_df %>%
  dplyr::select(datasetID, cluster1, OR)

OR_wide <- OR_selected %>%
  pivot_wider(names_from = datasetID, values_from = OR) %>%
  as.data.frame()

rownames(OR_wide) <- OR_wide$cluster1
OR_wide$cluster1 <- NULL

# Check how many cell types appear in at least one dataset
table(rowSums(!is.na(OR_wide)) >= 1)

# Visualize
library(ggplot2)
library(reshape2)
library(tibble)

# Log-transform OR for visualization
OR_wide_log <- log10(OR_wide)
OR_long <- OR_wide_log %>%
  tibble::rownames_to_column(var = "Cell_Type") %>%
  melt(id.vars = "Cell_Type", variable.name = "Dataset", value.name = "log10_OR")

# Classify OR direction for coloring
OR_long$Color_Group <- ifelse(is.na(OR_long$log10_OR), "NA",
                              ifelse(OR_long$log10_OR > 0, "High", "Low"))

# Reshape p-values for merging
pval_long <- pval_wide %>%
  tibble::rownames_to_column(var = "Cell_Type") %>%
  melt(id.vars = "Cell_Type", variable.name = "Dataset", value.name = "p_adj")

# Merge OR and p-values
heatmap_data <- left_join(OR_long, pval_long, by = c("Cell_Type", "Dataset"))

# Define significance (you may add a threshold)
heatmap_data$Significance <- ifelse(!is.na(heatmap_data$p_adj) & heatmap_data$p_adj < 0.05, "*", "")

# Count significant datasets per cell type
significant_count <- heatmap_data %>%
  group_by(Cell_Type) %>%
  summarise(
    Significant_Datasets = sum(Significance == "*", na.rm = TRUE),
    Present_Datasets = sum(!is.na(log10_OR))
  ) %>%
  mutate(Ratio = Significant_Datasets / Present_Datasets)

# Plot significance vs coverage with Ratio as point size
library(ggrepel)

ggplot(significant_count, aes(x = Present_Datasets, y = Significant_Datasets, 
                              size = Ratio)) +
  geom_point(color = "blue", alpha = 0.8) +
  geom_text_repel(aes(label = Cell_Type), size = 5, angle = 90, vjust = 0.5) +
  scale_size(range = c(2, 8)) +
  labs(title = "Significant Cell Type Enrichment Across Datasets",
       x = "Number of Datasets with Cell Type Present",
       y = "Number of Significant Enrichments (FDR < 0.05)",
       size = "Significance Ratio") +
  theme_minimal()
