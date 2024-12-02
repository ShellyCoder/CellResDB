# Clear the environment
rm(list = ls())
options(stringsAsFactors = F)

setwd("/sc/processed_data/")
file_list <- list.files(pattern = "rds$")

library(msigdbr) 
library(gplots)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(Seurat)
library(tools)

# Set up the output directory for enrichment analysis results
enrich_output_dir <- "enrich/"
dir.create(enrich_output_dir)     
dir.create("cluster_marker/")

# Loop through each file in the file list
for (file in file_list) {
  tryCatch({
    # Read the Seurat object
    data <- readRDS(file)
    
    # Create a directory for the current file's enrichment results
    op_dir <- paste0(enrich_output_dir, file_path_sans_ext(file))
    if (!dir.exists(op_dir)) {
      dir.create(op_dir)
    }
    
    # Get unique clusters
    clusters <- unique(data@meta.data$cluster)
    
    # Loop through each cluster
    for (cluster in clusters) {
      # Clean cluster name to remove non-alphanumeric characters
      cleaned_cluster <- gsub("[^[:alnum:][:space:]+]", " ", cluster)
      
      # Check if there are only two treatments and if there are enough samples
      test1 <- c("Pre", "Post")
      test <- unique(data@meta.data$treatment[data@meta.data$cluster == cluster])
      if (setequal(test, test1) & min(table(data@meta.data$treatment[data@meta.data$cluster == cluster])) > 3) {
        
        # Subset data for the current cluster
        cell_ids <- rownames(data@meta.data[data@meta.data[["cluster"]] == cluster, ])
        degs_data <- subset(data, cells = cell_ids)
        
        # Set the identity of the subset to treatment (Pre vs Post)
        Idents(degs_data) <- "treatment"
        
        # Perform differential expression analysis between Pre and Post treatments
        degs <- FindMarkers(degs_data, 
                            logfc.threshold = 0.25, 
                            min.pct = 0.1, 
                            only.pos = FALSE,
                            ident.1 = "Pre", ident.2 = "Post") %>%
          mutate(gene = rownames(.))  # Add gene names as a column
        
        # Filter the results for significant genes
        degs_fil <- degs %>% 
          filter(p_val < 0.05) %>%
          filter(abs(avg_log2FC) > 0.25)
        
        # Map gene symbols to ENTREZ IDs
        ids <- bitr(degs_fil$gene, 'SYMBOL', 'ENTREZID', 'org.Hs.eg.db')
        
        # Merge with ENTREZ IDs
        degs_fil <- merge(degs_fil, ids, by.x = 'gene', by.y = 'SYMBOL')
        degs_fil$cluster <- cluster
        
        # Save the results for each cluster
        write.csv(degs_fil, paste0("cluster_marker/", file_path_sans_ext(file), "-", cleaned_cluster, "_marker.csv"))
        
        # Sort the genes by log2 fold change
        degs_fil <- degs_fil[order(degs_fil$avg_log2FC, decreasing = TRUE), ]
        
        # Create a gene list based on log2FC values
        degs_list <- as.numeric(degs_fil$avg_log2FC)
        names(degs_list) <- degs_fil$ENTREZID
        
        # Select genes with large log2FC
        cluster_de <- names(degs_list)[abs(degs_list) > 1]
        
        # GO enrichment analysis
        cluster_ego <- enrichGO(cluster_de, OrgDb = "org.Hs.eg.db", ont = "BP", readable = TRUE, pvalueCutoff = 0.2)
        
        if (length(rownames(head(cluster_ego))) > 0) {
          # Save the GO results
          ego_bp <- data.frame(cluster_ego@result)
          write.csv(ego_bp, paste0(op_dir, "/", file_path_sans_ext(file), "_", cleaned_cluster, "_go.csv"))
          
          # Plot the GO bar chart
          p2 <- barplot(cluster_ego, font.size = 14, showCategory = 10, color = "pvalue") +
            theme(plot.margin = unit(c(1, 1, 1, 1), 'lines'))
          print(p2)
          
          # Save the plot
          ggsave(paste0(op_dir, "/", file_path_sans_ext(file), "_", cleaned_cluster, "_go.jpg"), plot = p2, width = 12, height = 10, dpi = 300, device = "jpeg", limitsize = FALSE)
        }
        
        # KEGG enrichment analysis
        cluster_ekg <- enrichKEGG(gene = cluster_de, organism = "hsa", pvalueCutoff = 0.2)
        
        if (length(rownames(head(cluster_ekg))) > 0) {
          # Save the KEGG results
          kegg_result <- data.frame(cluster_ekg@result)
          write.csv(kegg_result, paste0(op_dir, "/", file_path_sans_ext(file), "_", cleaned_cluster, "_kegg.csv"))
          
          # Plot the KEGG bar chart
          p3 <- barplot(cluster_ekg, font.size = 14, showCategory = 10, color = "pvalue") +
            theme(plot.margin = unit(c(1, 1, 1, 1), 'lines'))
          print(p3)
          
          # Save the plot
          ggsave(paste0(op_dir, "/", file_path_sans_ext(file), "_", cleaned_cluster, "_kegg.jpg"), plot = p3, width = 12, height = 10, dpi = 300, device = "jpeg", limitsize = FALSE)
        }
      }
    }
  }, error = function(e) {
    # Handle errors by printing the file and error message
    print(paste0(file, ": ", e))
  })
}



