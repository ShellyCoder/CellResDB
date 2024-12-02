# Clear the environment
rm(list = ls())
options(stringsAsFactors = F)

# Load necessary library
library(Seurat)

# Set working directory
setwd("/sc/data/")           
file_list <- list.files(pattern = "rds$")

# Create output directory if it doesn't exist
output_dir <- "/sc/processed_data/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Process each .rds file
sceList = lapply(file_list, function(x){
  print(x)
  
  # Load the Seurat object
  sce <- readRDS(x)
  
  # Calculate the percentage of mitochondrial gene expression
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
  
  # Filter cells based on feature count and mitochondrial percentage
  sce <- subset(sce, subset = nFeature_RNA > 200 & percent.mt < 20)
})

# Process each Seurat object and save the results
for (i in 1:length(sceList)) {
  # Normalize data
  sceList[[i]] <- NormalizeData(sceList[[i]], verbose = FALSE, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Find variable features
  sceList[[i]] <- FindVariableFeatures(sceList[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  
  # Scale the data for variable features
  sceList[[i]] <- ScaleData(sceList[[i]], features = VariableFeatures(sceList[[i]]))
  
  # Perform PCA
  sceList[[i]] <- RunPCA(sceList[[i]], features = VariableFeatures(object = sceList[[i]]))
  
  # Define output file name
  output_file <- paste0(output_dir, file_list[i])  
  
  # Save the processed Seurat object
  saveRDS(sceList[[i]], file = output_file)  
  print(paste("Saved:", output_file))  
}



