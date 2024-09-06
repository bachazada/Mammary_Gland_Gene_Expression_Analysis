# Mammary_Gland_Gene_Expression_Analysis
This repository contains analysis of gene expression profiles of basal stem-cell enriched cells (B) and committed luminal cells (L) in the mammary gland of virgin, pregnant, and lactating mice  includes differential expression, quality control, and visualization using heatmaps and multi-dimensional scaling (MDS) plots.

# Mammary Gland Gene Expression Analysis

This repository contains the scripts, data, and results for analyzing the gene expression profiles of basal stem-cell enriched cells (B) and committed luminal cells (L) in the mammary gland of virgin, pregnant, and lactating mice. The analysis includes quality control, filtering, normalization, and visualizations such as heatmaps and MDS plots.

## Files and Directories

- **R Scripts:**
  - `gene_expression_analysis.R`: Contains the code for loading data, filtering, normalization, and generating plots.
  
- **Data Files:**
  - `GSE60450_Lactation-GenewiseCounts.txt`: Gene-wise counts data for different cell types and statuses.
  - `SampleInfo.txt`: Information about the samples, including cell type and status.
  - `SampleInfo_Corrected.txt`: Corrected sample information.

- **Results:**
  - `bar_plot_library_sizes.pdf`: A bar plot showing the library sizes of each sample.
  - `distribution_log2_CPMs.pdf`: A box plot showing the distribution of log-transformed counts per million (CPM) values.
  - `MDS.pdf`: Multi-dimensional scaling (MDS) plots showing cell type and status differences.
  - `heatmap_most_variable.pdf`: A heatmap showing the samples' top 500 most variable genes.

- **Documentation:**
  - `Tutorial_4_NGS_analysis.pdf`: Contains step-by-step instructions for performing gene expression analysis using R.

## Usage

1. Clone the repository:
git clone https://github.com/yourusername/Mammary_Gland_Gene_Expression_Analysis.git


2. Open the R script `gene_expression_analysis.R` to run the analysis. Make sure to install the necessary R packages if not already installed:
```R
install.packages(c("edgeR", "limma", "Glimma", "gplots", "RColorBrewer", "org.Mm.eg.db", "GO.db"))
