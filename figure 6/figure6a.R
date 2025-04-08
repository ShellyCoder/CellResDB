# This script performs a 2×2 contingency table test (Chi-square or Fisher's exact test) 
# to assess the enrichment of specific cell types between response (R) and non-response (NR) 
# groups in post-treatment samples across multiple datasets and diseases.

rm(list = ls())
options(stringsAsFactors = FALSE)

# Load required libraries
library(readxl)
library(dplyr)

# Load cell type matching information
cellTable <- read_excel("data/revision_data/match.xlsx")
table(is.na(cellTable))

length(unique(cellTable$new))
length(unique(cellTable$old))

# Load all meta data files
files <- list.files("data/revision_data/meta_data/", pattern = "\\.txt$", full.names = TRUE)
df_list <- lapply(files, read.table, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Extract dataset names from file names
dataNames <- unlist(purrr::map(files, basename)) %>% stringr::str_remove("_meta_data.txt$")

# Combine metadata into a unified data frame
df_list <- lapply(1:length(df_list), function(i) {
  tmp_df <- df_list[[i]]
  tmp_df <- tmp_df[, c("patient", "treatment", "Response", "cluster1")]
  tmp_df$datasetID <- dataNames[i]
  return(tmp_df)
})

df_combined <- do.call(rbind, df_list)

# Add tissue and disease annotations
datasetInfo <- read.table("data/revision_data/tissue_Search.txt", sep = "\t", 
                          stringsAsFactors = FALSE, header = TRUE)
df_combined <- dplyr::left_join(df_combined, datasetInfo, by = "datasetID")

# Subset to post-treatment samples
postCell_df <- df_combined[df_combined$treatment == "Post", ]

# Sort by tissue and disease
postCell_sorted <- postCell_df %>%
  arrange(tissue, standardDisease)

diseases <- unique(postCell_sorted$standardDisease)
final_results <- list()

for (d in diseases) {
  
  print(d)
  postCell_tmp <- postCell_sorted[postCell_sorted$standardDisease == d, ]
  dataset_tmp <- unique(postCell_tmp$datasetID)
  
  # Only include datasets where both response types (R and NR) are present
  results_list <- lapply(dataset_tmp, function(x) {
    print(x)
    tmp_df <- postCell_tmp[postCell_tmp$datasetID == x, ]
    
    if (length(unique(tmp_df$Response)) == 2) {
      
      tmp_df$Response <- factor(tmp_df$Response, levels = c("NR", "R"))
      
      # Perform 2×2 contingency table test for each cell type
      cell_type_counts <- table(tmp_df$cluster1, tmp_df$Response)
      p_values <- c()
      OR_values <- c()
      
      # Filter out low count cell types
      target_cellType <- rownames(cell_type_counts)[rowSums(cell_type_counts < 20) == 0]
      if (length(target_cellType) < 1) return(NULL)
      
      for (cell_type in target_cellType) {
        
        # Construct 2×2 contingency table
        cell_counts <- tmp_df %>%
          mutate(is_cell = ifelse(cluster1 %in% cell_type, "Cell", "Non Cell")) %>%
          count(is_cell, Response) %>%
          tidyr::spread(Response, n, fill = 0)
        
        if (all(c("Cell", "Non Cell") %in% cell_counts$is_cell)) {
          contingency_table <- as.matrix(cell_counts[, 2:3])
          rownames(contingency_table) <- cell_counts$is_cell
        } else {
          stop("Error in 2×2 contingency table construction.")
        }
        
        # Perform statistical test
        if (any(contingency_table < 5)) {
          test_result <- fisher.test(contingency_table)
          OR_tmp <- test_result$estimate
        } else {
          test_result <- chisq.test(contingency_table)
          OR_tmp <- fisher.test(contingency_table)$estimate
        }
        
        # Store p-values and odds ratios (OR > 1 means enriched in NR group)
        p_values[cell_type] <- test_result$p.value
        OR_values[cell_type] <- ifelse(is.null(OR_tmp), NA, OR_tmp)
      }
      
      # Adjust p-values using FDR correction
      p_adj <- p.adjust(p_values, method = "fdr")
      
      result_df <- data.frame(
        datasetID = x,
        cluster1 = names(p_values),
        p_value = p_values,
        p_adj = p_adj,
        OR = OR_values
      )
      return(result_df)
      
    } else {
      return(NULL)
    }
  })
  
  names(results_list) <- dataset_tmp
  final_results <- c(final_results, results_list)
}

# Clean and merge results
clean_result <- final_results[!sapply(final_results, is.null)]
related_df <- do.call("rbind", clean_result)
