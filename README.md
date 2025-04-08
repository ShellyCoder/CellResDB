# CellResDB
CellResDB: deciphering cancer therapy resistance via patient-level single-cell transcriptomics

This repository contains scripts used for data collection and data analysis in this work.

## Annotation
- **`1.Data Preprocessing/`**:
  - **`01_Data preprocess.R`**: This script is responsible for preprocessing single-cell transcriptomic data. It includes steps such as normalization, filtering, and quality control to prepare the data for further analysis.

- **`2.Downstream Analysis/`**:
  - **`02_UMAP and DEGs.R`**: This script performs UMAP analysis for dimensionality reduction and visualizes the results. It also conducts differential expression analysis to identify genes that are significantly expressed in different cell types or conditions.
  - **`03_TME composition.R`**: This script analyzes the tumor microenvironment (TME) composition, providing insights into the cellular makeup and interactions within the tumor ecosystem.
  - **`04_Functional enrichment.R`**: This script conducts functional enrichment analysis to identify biological processes and pathways that are affected in the differentially expressed genes.
  - **`05_Cell communication.R`**: This script analyzes intercellular communication, exploring how different cell types interact with each other within the tumor microenvironment.
  - **`06_GSEA.R`**: This script performs GSEA enrichment analysis on paired response (R) and non-response (NR) samples.
- **`figure3/`**: Scripts for analysis in Figure 3.
  - **`figure3a.R`**: Creates a dendrogram with heatmap annotations showing sample attributes for cancer subtypes and tissues.
  - **`figure3b.R`**: Constructs a star tree visualization of sample relationships, organized by tissue and response types.
  - **`figure3c.R`**: Produces a bar plot to show the frequency distribution of treatment types sorted alphabetically.
  - **`figure3d.R`**: Visualizes individual drug usage distribution across therapeutic regimens with a horizontal bar chart.
  - **`figure3e.R`**: Summarizes and visualizes the number of unique cell types across datasets using a binned bar plot.
  - **`figure3f.R`**: Plots the distribution of cell counts (log10-transformed) using binned bar charts.

- **`figure6/`**: Scripts for analysis in Figure 6.
  - **`figure6a_and_6b.R`**: Analyzes relative proportions of treatment responses (NR and R) and constructs hierarchical clustering trees for cancer types.
  - **`figure6c.R`**: Generates a bubble plot to visualize the relationship between cancer types and treatment combinations.

## Usage
Clone this repository:
   ```bash
   git clone https://github.com/ShellyCoder/CellResDB.git
