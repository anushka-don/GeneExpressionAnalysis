library(dbplyr)
library(tidyverse)
library(gridExtra)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(ggrepel)
library(DESeq2)
library(DT)
library(igraph)
library(matrixStats)

# Function for reading Series Matrix file for the given dataset
read_series_mat <- function(file) {
  data <- read.table(text = grep("!Sample", readLines(file), value = TRUE), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  # Remove the prefix "!Sample_" from each title
  row_names <- gsub("!Sample_", "", data[[1]])
  row.names(data) <- make.unique(row_names)
  
  # Remove the first column (which is now used as row names)
  data <- data[-1]
  
  # Set column names from the second row
  col_names <- data[2, ]
  colnames(data) <- col_names
  return(data)
}


# Function for summarizing sample information of the dataset
modify_series_mat <- function(data) {
  dataf <- t(data)

  # Select columns where names start with "characteristics_ch1"
  selected_columns <- grep("^characteristics_ch1", colnames(dataf), value = TRUE)
  
  # Subset the data to include only those columns
  dataf <- dataf[, selected_columns, drop = FALSE]
  
  # Extract column names from the first row
  col_names <- sapply(dataf[1,], function(x) gsub(":.*", "", x))

  # Remove the first row and set column names
  colnames(dataf) <- col_names
  
  return(dataf)
}


# Function for making the summary table 
summarize_series_mat <- function(dataf) {
  # Count the number of samples
  total_samples <- nrow(dataf)
  
  # Make Characteristic and Value columns for the summary of the dataset
  summary_table <- dataf %>%
    as.data.frame() %>%
    pivot_longer(cols = everything(), names_to = "Characteristic", values_to = "Value") %>%
    mutate(Value = gsub(".*: ", "", Value)) %>%
    group_by(Characteristic, Value) %>%
    summarise(Count = n(), .groups = "drop")
  
  # Format the summary table
  summary_table <- summary_table %>%
    group_by(Characteristic) %>%
    summarise(Value = paste(paste0(Value, " (", Count, ")"), collapse = ", ")) %>%
    ungroup() %>% 
    bind_rows(
      data.frame(Characteristic = "Total Samples", Value = as.character(total_samples)),.) 
  return(summary_table)
}


# Function for sample information
sample_info <- function(data,smpl_name) {
  
  # Remove the second row as it is now used as column names
  data <- data[-2, ]
  
  return(data.frame('SampleCharacteristics' = row.names(data), Value = data[[smpl_name]]))
}


# Function to create pie chart for a given column
create_pie_chart <- function(data, column_name) {
  # Count occurrences
  count_data <- data %>%
    as.data.frame() %>%
    pivot_longer(cols = everything(), names_to = "Characteristic", values_to = "Value") %>%
    mutate(Value = gsub(".*: ", "", Value)) %>%
    group_by(Characteristic, Value) %>%
    summarise(Count = n(), .groups = "drop")
  
  filtered_data <- count_data %>% filter(Characteristic == column_name)
  
  # Create pie chart
  plt <- ggplot(filtered_data, aes(x = "", y = Count, fill = Value)) +
    geom_bar(stat = "identity") +
    coord_polar(theta = "y") +
    labs(title = paste("Pie Chart of", column_name), x = NULL, y = NULL) +
    theme_void() +
    scale_fill_brewer(palette = "Set1") # Optional color palette
  
  return(plt)
}

# Function to combine all the gene files
combine_gene_files <- function(files) {
  # Read the first file to initialize the result dataframe
  f1 <- read.delim(files$datapath[1], sep = '\t', header = FALSE)
  counts_df <- data.frame(row.names = f1[[1]])  # Use the first column as row names
  
  for (i in seq_len(nrow(files))) {  # Iterate over each uploaded file
    # Read the current file
    data <- read.delim(files$datapath[i], sep = '\t', header = FALSE)
    
    # Extract the third column
    third_column <- data[[3]]  # Adjust index if necessary
    
    # Extract the file name and create a column name
    original_name <- basename(files$name[i])  # Use the 'name' column for the original file name
    column_name <- sub("_.*", "", original_name)  # Removes everything after and including the first underscore
    
    # Add it to counts_df with the modified filename as the column name
    counts_df[[column_name]] <- third_column
  }
  
  counts_df <- rownames_to_column(counts_df,var = 'gene')
  return(counts_df)
}

# Normalize raw count data to counts per million WITH pseudocounts using the following formula:
#'     count / (sample_library_size/10^6)
get_library_size <- function(count_data) {
  # Calculate total read counts for each sample, by adding all the values in a column
  total_counts <- count_data %>%
    summarise(across(-gene, sum, na.rm = TRUE)) 
  
  # visualise the values by making the tibble longer
  result <- total_counts %>%
    pivot_longer(cols = everything(), names_to = "sample", values_to = "value")
  
  return(result)
}
normalize_by_cpm <- function(count_data) {
  # Make the output from the library_size wider (horizontal)
  library_sum <- get_library_size(count_data) %>%
    pivot_wider(names_from = 'sample', values_from = 'value')
  
  # Convert the library_sum to a numeric vector
  library_sizes <- as.numeric(library_sum[1, ])
  
  # Normalize counts by dividing each column by the respective library size in millions
  norm_counts <- count_data %>% 
    column_to_rownames(var = 'gene') %>% 
    sweep(2, library_sizes * 1e-6, FUN = "/") %>%
    rownames_to_column(var = 'gene')
  
  norm_counts <- as_tibble(norm_counts)
  
  return(norm_counts)
}

# variance filter for the genes
filter_var_genes <- function(verse_counts, variance) {
  # Filter the genes whose variance is above given value, using filter
  filtered <- verse_counts %>%
    rowwise() %>%
    filter(var(c_across(-gene)) > variance) %>%
    ungroup()
  
  return(filtered)
}

# Filtering based on non-zero genes for samples
filter_non_zero_genes <- function(counts_df, threshold) {
  filtered_df <- counts_df %>%
    rowwise() %>%
    mutate(non_zero_count = sum(c_across(starts_with("GSM")) > 0)) %>%
    filter(non_zero_count >= threshold) %>%
    select(-non_zero_count)
  
  return(filtered_df)
}

# Function for Summary of Counts file
summary_counts <- function(counts_df, filtered_df) {
  # Get dimensions of the original and filtered dataframes
  genes_all <- dim(counts_df)[1]  # Number of rows (genes) in counts_df
  samples_all <- dim(counts_df)[2] # Number of columns (samples) in counts_df
  
  genes_fil <- dim(filtered_df)[1]  # Number of rows (genes) in filtered_df
  samples_fil <- dim(filtered_df)[2] # Number of columns (samples) in filtered_df
  
  # Calculate the number of genes passing the filter and not passing the filter
  genes_passing_filter <- genes_fil
  genes_not_passing_filter <- genes_all - genes_passing_filter
  
  # Calculate the percentage of genes passing and not passing the filter
  perc_passing_filter <- (genes_passing_filter / genes_all) * 100
  perc_not_passing_filter <- (genes_not_passing_filter / genes_all) * 100
  
  # Create a summary dataframe
  summary_df <- data.frame(
    Metric = c("Total Genes", "Total Samples", "Genes Passing Filter", 
               "Genes Passing Percentage (%)", "Genes Not Passing Filter", 
               "Genes Not Passing Percentage (%)"),
    Value = c(genes_all, samples_all - 1, genes_passing_filter, 
              perc_passing_filter, genes_not_passing_filter, 
              perc_not_passing_filter)
  )
  
  return(summary_df)
}


# Function for scatter plots 
calculate_metrics <- function(counts_df) {
  metrics_df <- counts_df %>%
    rowwise() %>%
    mutate(
      Median_Count = median(c_across(starts_with("GSM")), na.rm = TRUE),
      Variance_Count = var(c_across(starts_with("GSM")), na.rm = TRUE),
      Number_of_Zeros = sum(c_across(starts_with("GSM")) == 0)
    ) %>%
    ungroup()
  
  return(metrics_df)
}


# Function for Differential Expression of counts file
run_deseq <- function(counts_df, metadata_table, design_formula) {
  
  # Ensure 'gene' column is used as rownames in counts_df
  counts_df <- column_to_rownames(counts_df, var = 'gene')
  
  metadata_table <- metadata_table[match(colnames(counts_df), rownames(metadata_table)), ]
  
  # Ensure counts are integers
  counts_df <- round(counts_df)
  
  print(ncol(counts_df))
  print(nrow(metadata_table))
  
  # Create DESeq2 dataset using dynamic design formula
  dds <- DESeqDataSetFromMatrix(
    countData = counts_df,
    colData = metadata_table,
    design = as.formula(design_formula) # Use the user-specified design
  )
  
  # Filter out rows with low counts
  dds <- dds[rowSums(counts(dds)) > 1, ]
  
  # Run DESeq2 to normalize and perform differential expression analysis
  dds <- DESeq(dds)
  
  # Extract results
  res <- results(dds) %>%
    as.data.frame() %>%
    rownames_to_column(var = 'gene')
  
  return(res)
}

# Function for labeling the results
label_res <- function(deseq2_res, padj_threshold) {
  deseq2_res <- deseq2_res %>%
    mutate(volc_plot_status = case_when(
      padj < padj_threshold & log2FoldChange > 0 ~ "UP", 
      padj < padj_threshold & log2FoldChange < 0 ~ "DOWN",
      TRUE ~ "NS"
    ), .after = gene) %>%
    as_tibble()
  
  return(deseq2_res)
}


# Function for plotting DESeq results - pvalue
plot_pvals <- function(labeled_results) {
  p <- ggplot(labeled_results, aes(x = pvalue)) +
    geom_histogram(binwidth = 0.02, position = 'identity',fill = "lightblue", color = "black") +
    labs(title = "Histogram of Unadjusted P-Values", x = "pvalue", y = "count") +
    theme_minimal()
  
  return(p)
}


# Function for plotting DESeq results - Log2FoldChange
plot_log2fc <- function(labeled_results, padj_threshold) {
  # Filter for significant genes based on padj threshold
  significant_genes <- labeled_results %>%
    filter(padj < padj_threshold)
  
  # Plot a histogram of log2 fold changes for significant genes
  p <- ggplot(significant_genes, aes(x = log2FoldChange)) +
    geom_histogram(binwidth = 0.2, fill = "lightgreen", color = "black") + 
    labs(title = paste("Histogram of Log2 Fold Changes (Significant at padj <", padj_threshold, ")"), 
         x = "Log2 Fold Change", 
         y = "Frequency") +
    theme_minimal()
  
  return(p)
}


# Function for plotting Volcano plot - 
plot_volcano <- function(labeled_results) {
  p <- ggplot(labeled_results, aes(x = log2FoldChange, y = -log10(padj), color = volc_plot_status)) +
    geom_point(size = 2) +  # Adjust alpha and size for visibility
    labs(title = "Volcano Plot of DESeq2 Results",
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-value") +
    theme_minimal() +
    scale_color_manual(values = c("UP" = "blue", "DOWN" = "red", "NS" = "green")) +  # Color settings
    theme(legend.title = element_blank())
  
  return(p)
}
