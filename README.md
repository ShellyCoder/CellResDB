# CellResDB
CellResDB: deciphering cancer therapy resistance via patient-level single-cell transcriptomics

This repository contains scripts used for data collection and data analysis in this work.

## Directory Structure
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
