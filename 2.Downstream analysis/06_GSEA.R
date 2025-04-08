# This script performs GSEA enrichment analysis on paired response (R) and non-response (NR) samples
# within the same dataset for a specific cell type (Monocyte/Macrophage).
# The goal is to identify treatment-associated transcriptional signatures using hallmark gene sets.

rm(list = ls())
options(stringsAsFactors = FALSE)

# Load necessary libraries
library(readxl)
library(dplyr)
library(clusterProfiler)
library(enrichplot)

# Path to MSigDB hallmark gene set (download required in GMT format)
gmt_file <- "./h.all.v2024.1.Hs.symbols.gmt"
gene_sets <- read.gmt(gmt_file)  # Load gene sets

# List of dataset identifiers
marcophagy_files <- c("GSE120575_melanoma_2", "GSE120575_melanoma_1", "GSE123813_BCC_1", 
                      "GSE206325_HCC_2", "GSE169246_TNBC_2", "GSE161801_MM_2", 
                      "GSE199333_AML", "GSE111014_CLL")

# Initialize list to store GSEA results
GSEA_list <- list()

for (f in marcophagy_files) {
  
  # Load Seurat object
  tmp_SeuratObj <- readRDS(paste0("./my_directory/", f, ".rds"))
  
  # Subset to post-treatment samples and Monocyte/Macrophage cluster
  post_cells <- subset(tmp_SeuratObj, subset = treatment == "Post")
  post_cells <- subset(post_cells, subset = cluster1 == "Monocyte/Macrophage")
  post_cells@meta.data$Response <- factor(post_cells@meta.data$Response)
  
  # Set identity to Response status (R vs NR)
  Seurat::Idents(post_cells) <- post_cells@meta.data$Response
  print(table(post_cells@meta.data$Response))
  
  # Differential expression analysis between R and NR
  deg_df <- Seurat::FindMarkers(post_cells, ident.1 = "R", ident.2 = "NR", 
                                        only.pos = FALSE, logfc.threshold = 0, 
                                        min.pct = 0.1, test.use = "wilcox")
  
  # Filter and sort by log2FC
  deg_filtered <- deg_df %>%
    filter(!is.na(avg_log2FC)) %>%
    arrange(desc(avg_log2FC))
  
  # Clear memory
  rm(tmp_SeuratObj, post_cells, deg_df)
  gc()
  
  # Prepare gene list for GSEA
  gene_list <- deg_filtered$avg_log2FC
  names(gene_list) <- rownames(deg_filtered)
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # Perform GSEA using hallmark gene sets
  gsea_results <- GSEA(geneList = gene_list, TERM2GENE = gene_sets, 
                       pvalueCutoff = 0.05)
  
  # Store results if available
  if (!is.null(gsea_results)) {
    GSEA_list[[f]] <- gsea_results@result
  }
}
