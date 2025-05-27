Website - https://gene-expression-analysis-geo.shinyapps.io/final_project_r/
GEO Dataset - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266358  
Research article Link - https://www.nature.com/articles/s41419-024-07152-0#Abs1

# README

## Overview

This Shiny application, titled **"Studying Genes involved in Alzheimer's Disease Progression,"** is a tool designed for researchers and analysts working on Alzheimer's disease. It provides an interactive platform to analyze gene expression data, study sample metadata, perform primary analysis, and conduct differential expression and network analysis.

The application is divided into multiple tabs, each catering to specific aspects of gene expression and sample analysis. Users can upload data files, visualize results through plots, and explore complex relationships between genes through network analyses.

---

## Features

### 1. **Samples Tab**
This tab focuses on metadata and sample-related information.
- **File Upload:** Users can upload a Series Matrix File in `.csv`, `.tsv`, or `.txt` format.
- **Tabs within Samples Panel:**
  - **Summary:** Displays a summary of the uploaded sample data.
  - **Sample Info:** A detailed, interactive table showing sample metadata.
  - **Plots:** Visualizations such as pie charts of sample characteristics.
- **Error Handling:** If no file is uploaded, an appropriate error message is displayed.

### 2. **Counts Tab**
This tab is used for primary analysis of gene expression data.
- **File Upload:**
  - Upload one or more gene count files in `.csv`, `.tsv`, or `.txt` format.
- **Analysis Options:**
  - Filter genes by variance percentile and minimum number of non-zero samples.
  - Download normalized counts for further analysis.
- **Tabs within Counts Panel:**
  - **Filter Summary:** Displays a summary of filtered genes.
  - **Scatter Plots:** Shows scatter plots of variance and non-zero samples.
  - **Heatmap:** A clustered heatmap for gene and sample visualization.
  - **PCA Analysis:** Perform PCA on the data and visualize principal components with scatter plots and beeswarm plots.
    - Customize PCA plots by selecting principal components for the x-axis and y-axis.
    - Specify the number of top components for additional visualization.

### 3. **DE (Differential Expression Analysis) Tab**
This tab performs differential expression analysis using DESeq2.
- **Analysis Options:**
  - Select a design factor (e.g., `~ treatment`, `~ genotype`, or combined factors).
  - Adjust the p-value threshold for significance.
- **Tabs within DE Panel:**
  - **DESeq Results:** Interactive table displaying differential expression results.
  - **Plots:**
    - Histogram of raw p-values.
    - Log2 fold change histogram.
    - Volcano plot for visualizing differential expression.

### 4. **Network Tab**
This tab allows users to explore gene networks and correlations.
- **File Upload:**
  - Upload a normalized counts matrix.
- **Analysis Options:**
  - Enter a list of genes for network analysis.
  - Set a correlation threshold to define network edges.
  - Explore shortest paths between selected vertices.
- **Tabs within Network Panel:**
  - **Clustered Heatmap:** Visualize gene clusters.
  - **Correlation Network:** A graphical representation of gene correlations.
  - **Network Metrics:** Summary of network parameters.
  - **Shortest Path:** Compute and display the shortest path between selected vertices.

---

## Getting Started

### Prerequisites
- Install the `shiny` package in R.
- Include additional dependencies as required by `functions.R`.

### Installation
1. Clone or download this repository.
2. Ensure the file `functions.R` is in the same directory as the Shiny app file.
3. Launch the app by running the following command in R:
   ```R
   shiny::runApp()
   ```

### File Formats Supported
- `.csv`
- `.tsv`
- `.txt`

### Configurations
- **Maximum File Size:** The app supports uploads up to 30 MB. This can be adjusted using the `shiny.maxRequestSize` option.

---

## Usage
1. Navigate through the tabs to perform different analyses.
2. Upload the required data files in the appropriate sections.
3. Use the interactive sliders, buttons, and inputs to customize the analysis.
4. Visualize and download results as needed.

---

## Limitations
- Ensure that the input files are formatted correctly to avoid errors.
- Proper annotations and metadata are required for meaningful analysis.

---

## To-Do
The application is still under development. Future updates may include:
- Enhanced error handling.
- Additional plotting, network creation and visualization options.
- Integration with external APIs for data retrieval.

---

## Contact
For questions or issues, please contact the developer at [d13anushka@gmail.com].
