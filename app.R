library(shiny)
source('functions.R')

# Increase the maximum file upload size
options(shiny.maxRequestSize = 30 * 1024^2)


# UI for the Application
ui <- fluidPage(

    # Application title
    titlePanel("Studying Genes involved in Alzheimer's Disease Progression"),
    
    ## Tabs
    tabsetPanel(
      
      # Tab 1 - Samples containing summary of the samples metadata, information about the study samples and Pie charts 
      tabPanel("Samples",
               sidebarLayout(
                 sidebarPanel(
                   # File input button
                   fileInput(inputId = 'seriesmat',
                             label = 'Upload your Series Matrix File here',
                             accept = c('.csv','.tsv','.txt')),
                   
                   actionButton(inputId = "submit", label = "Submit")),
                 
                 mainPanel(
                   tabsetPanel(
                     # First Tab - Summary of the samples uploaded
                     tabPanel("Summary",
                       tableOutput(outputId = 'summary'),
                           # In case no file is uploaded
                       textOutput(outputId = 'error_message_series1')),
                     
                     # Second Tab - Information about the sample
                     tabPanel('Sample Info',
                              DTOutput(outputId = 'sample_information_table')),
                     
                     # Third Tab - Pie charts of the characteristics of the samples
                     tabPanel('Plots',
                              plotOutput(outputId = 'plot_pie'))
                     )
                   )
                 )
               ),
      
      # Tab 2 - Primary analysis of the Gene expression data 
      tabPanel("Counts",
               sidebarLayout(
                 sidebarPanel(
                   
                   # Upload all the gene counts files
                   fileInput(inputId = 'inputfiles',
                             label = 'Upload your Gene Counts File(s) here',
                             multiple = TRUE,
                             accept = c('.csv', '.tsv', '.txt')),
                   
                   sliderInput(inputId = 'variance',
                               label = 'Variance percentile',
                               min = 0,
                               max = 100,
                               value = 25),
                   
                   sliderInput(inputId = 'non_zero_samples',
                               label = 'Minimum number of non-zero samples',
                               min = 0,
                               max = 32,
                               value = 5),
                   
                   actionButton(inputId = "submit2", label = "Submit"),
                   
                   # Download normalized counts, in case you don't have the data
                   downloadButton("download_norm_counts", "Download Normalized Counts")
                   
                 ),
                 
                 mainPanel(
                   tabsetPanel(
                     
                     # First tab - Summary of all the filtered genes
                     tabPanel("Filter Summary",
                              tableOutput(outputId = 'summary_genes'),
                              # In case no file(s) are uploaded
                              textOutput(outputId = 'error_message_series2')
                     ),
                     
                     # Second Tab - Scatter plots
                     tabPanel('Scatter plots',
                              plotOutput(outputId = 'scatter_plot_var'),
                              plotOutput(outputId = 'scatter_plot_non_zero'),
                     ),
                     
                     # Third Tab - Heatmap - clustered for genes and samples visualization
                     tabPanel('Heatmap',
                              plotOutput(outputId = 'heat_map')
                     ),
                     
                     # Fourth Tab - Performing PCA and visualization of PCs
                     tabPanel(
                       "PCA",
                       sidebarLayout(
                         sidebarPanel(
                           h4("PCA Options"),
                           selectInput(
                             inputId = "pc_x",
                             label = "Select X-axis Principal Component:",
                             choices = NULL,
                             selected = NULL
                           ),
                           selectInput(
                             inputId = "pc_y",
                             label = "Select Y-axis Principal Component:",
                             choices = NULL,
                             selected = NULL
                           ),
                           numericInput(
                             inputId = "top_n",
                             label = "Number of Top Principal Components for Beeswarm Plot:",
                             value = 5,
                             min = 1
                           ),
                           actionButton(
                             inputId = "plot_beeswarm",
                             label = "Plot PCA box plot"
                           )
                         ),
                         mainPanel(
                           plotOutput(outputId = "pca_scatter"),
                           plotOutput(outputId = "pca_beeswarm")
                         )
                       )
                     )
                   )
                 )
               )
      ),
      
      
      # Tab 3 - Performing Differential Expression Analysis using DESeq2
      tabPanel("DE",
               
               sidebarLayout(
                 sidebarPanel(
                   
                   selectInput("design_factor", "Select Design Factor", 
                               choices = c("~ treatment", "~ genotype", "~ treatment + genotype")),
                   
                   sliderInput(input = 'padj_threshold',
                               label = 'Please Select the Adjusted P value',
                               min = 0,
                               max = 1,
                               value = 0.1),
                   
                   actionButton("run_analysis", "Run DESeq")
                 ),
                 
                 mainPanel(
                   tabsetPanel(
                     
                     # First Tab - DESeq Results
                     tabPanel(
                       title = 'DESeq Results',
                       DTOutput(outputId = 'deseq_results'),
                     ),
                     
                     # Second Tab - Plots in different tabs
                     tabPanel(
                       title = 'Plots',
                       tabsetPanel(
                         tabPanel(title = 'Raw p-values',
                                  plotOutput(outputId = 'pvalue_hist'
                                  )),
                         tabPanel(title = 'Log2FoldChanges',
                                  plotOutput(outputId = 'log2fc_hist'
                                  )),
                         tabPanel(title = 'Volcano plot',
                                  plotOutput(outputId = 'volcano_scatter'
                                  ))
                       )
                     )
                   )
                  )
               )
               ),
      
      # Tab 4 - Network
      tabPanel("Network",
               sidebarLayout(
                 sidebarPanel(
                   
                   fileInput("datafile", "Upload Normalized Counts Matrix:", accept = "c(.txt, .csv, .tsv)"),
                   
                   textAreaInput(
                     "gene_list",
                     "Enter Gene Names (one per line):",
                     rows = 10,
                     placeholder = "Enter gene names here..."),
                   
                   sliderInput(
                     "cor_threshold",
                     "Minimum Correlation Threshold:",
                     min = 0, max = 1, value = 0.4, step = 0.01),
                   
                   selectizeInput("vertex1", "Select Vertex 1:", choices = NULL),
                   selectizeInput("vertex2", "Select Vertex 2:", choices = NULL),
                   
                   actionButton("shortest_path_btn", "Show Shortest Path")),
                 
                 mainPanel(
                   
                   tabsetPanel(
                     # First Tab - Clustered Heatmap
                     tabPanel("Clustered Heatmap", 
                              plotOutput("heatmap")),
                     
                     # Second Tab - Correlation network of the selected genes
                     tabPanel("Correlation Network", 
                              plotOutput("network_plot")),
                     
                     # Third tab - Network metrics - Parameters used to create the network
                     tabPanel(
                       "Network Metrics",
                       dataTableOutput("network_metrics")),
                     
                     # Fourth Tab - Computing the shortest path
                     tabPanel("Shortest Path", 
                              verbatimTextOutput("shortest_path"))
                   )
                 )
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  
  ## First Tab output - Summary and Information of the samples 
  # Defining metadata table as a reactive variable, to be used in DEseq tab as well
  metadata_table <- reactiveVal(NULL)
  
  
  observeEvent(input$submit, {
    # Check if files are uploaded, return error otherwise
    if (is.null(input$seriesmat)) {
      output$error_message_series1 <- renderText("Error: Please upload the Series Matrix file before submitting.")
      output$summary <- renderTable(NULL)
      output$dynamic_select <- renderUI(NULL)
      output$sample_information_table <- renderTable(NULL)
    } 
    
    
    tryCatch({
      
      # Clear any previous error messages
      output$error_message_series1 <- renderText("")
      series_mat <- modify_series_mat(read_series_mat(input$seriesmat$datapath))
      
      metadata_table(series_mat)
      
      # First Sub Tab - Samples summary
      output$summary <- renderTable({
        summarize_series_mat(modify_series_mat(read_series_mat(input$seriesmat$datapath)))
      })
      

      # Second sub Tab - Sortable Sample Information as a DataTable
      output$sample_information_table <- renderDT({
        
        datatable(metadata_table(),
                  options = list(
                  pageLength = 5,        
                  lengthMenu = c(5, 10, 20), 
                  autoWidth = TRUE,      
                  dom = 'Bfrtip',
                  buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
                  ),
                  rownames = TRUE,
                  filter = 'top'
          )
      })
      
      # Third Sub Tab - Visualization of sample data as a Venn Diagram
      output$plot_pie <- renderPlot({
        req(input$seriesmat)
        
        # Store the modified series matrix for ease of reading for pie charts
        mod_reads <- modify_series_mat(read_series_mat(input$seriesmat$datapath))
        
        # Make Pie chart for all the characteristics to understand the distribution of characteristics in samples
        pie_charts <- lapply(colnames(mod_reads), function(col) {
          create_pie_chart(mod_reads, col)
        })
        
        grid.arrange(grobs = pie_charts, ncol = 2)
      })
      
    }, error = function(e) {
      output$error_message_series1 <- renderText("Oops! Wrong file. Please check the file.")
    })
  })
  
  
  
  ## Second Tab output - Counts data visualizer
  counts_file <- reactiveVal(NULL)
  filtered_genes <- reactiveVal(NULL)
  observeEvent(input$submit2, {
    # Check if files are uploaded, otherwise return an error
    if (is.null(input$inputfiles)) {
      output$error_message_series2 <- renderText("Error: Please upload the gene count files before submitting.")
      output$summary_genes <- renderTable(NULL)
      output$scatter_plots <- renderPlot(NULL)
      output$heat_map <- renderPlot(NULL)
      output$pca_plot <- renderPlot(NULL)
    }
    
    tryCatch({
      output$error_message_series2 <- renderText("")
      
      # Combine the uploaded counts files for all the samples and normalize the counts by counts per million method
      counts <- normalize_by_cpm(combine_gene_files(input$inputfiles))
      
      # Save counts file as a reactive value
      counts_file(counts)
      
      #filter the genes based on the criteria by the user
      filtered_genes(filter_non_zero_genes(filter_var_genes(counts_file(), input$variance),input$non_zero_samples))
      
      metrics_all <- calculate_metrics(counts_file())
      metrics_filtered <- calculate_metrics(filtered_genes())
      
      # Add a column to indicate if the gene is filtered or not
      metrics_all <- metrics_all %>%
        mutate(Filter_Status = "Not Passing")
      metrics_filtered <- metrics_filtered %>%
        mutate(Filter_Status = "Passing")
      combined_metrics <- bind_rows(metrics_all, metrics_filtered)

      output$download_norm_counts <- downloadHandler(
        filename = function() {
          paste("norm_counts.csv", sep = "")  # Name of the downloaded file
        },
        content = function(file) {
          write.csv(counts_file(), file, row.names = FALSE)  # Write the dataset to the file
        }
      )
      
      # Sub Tab 1 - Summary of the genes after applying the filters
      output$summary_genes <- renderTable({
        
        summary_counts(counts_file(), filtered_genes())
      })
      
      # Sub Tab 2 - Scatter plots
      output$scatter_plot_var <- renderPlot({
        ggplot(combined_metrics, aes(x = log10(Median_Count + 1), y = log10(Variance_Count + 1), color = Filter_Status)) +
          geom_point(alpha = 0.7) +
          scale_color_manual(values = c("Passing" = "darkblue", "Not Passing" = "lightgray")) +
          theme_minimal() +
          labs(
            title = "Median Count vs Variance",
            x = "Log10(Median + 1)",
            y = "Log10(Variance + 1)",
            color = "Filter Status"
          )
        })
      
      output$scatter_plot_non_zero <- renderPlot({
        # Scatter Plot: Median Count vs Number of Zeros
        ggplot(combined_metrics, aes(x = log10(Median_Count + 1), y = Number_of_Zeros, color = Filter_Status)) +
          geom_point(alpha = 0.7) +
          scale_color_manual(values = c("Passing" = "brown", "Not Passing" = "orange")) +
          theme_minimal() +
          labs(
            title = "Median Count vs Number of Zeros",
            x = "Log10(Median + 1)",
            y = "Number of Zeros",
            color = "Filter Status"
          )
        })
      
      
      # Sub Tab 3 - Heat Map
      output$heat_map <- renderPlot({
        matrix_data <- filtered_genes() %>%
          column_to_rownames(var = 'gene')
        
        matrix_data <- log10(matrix_data + 1)
        
        # Generate heatmap
        pheatmap(
          matrix_data,
          color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
          cluster_rows = TRUE,
          cluster_cols = TRUE,
          main = "Clustered Heatmap of Filtered Counts",
          legend = TRUE,
          legend_title = "Log10(Expression + 1)",
          show_rownames = FALSE
        )
      })
      
      # Sub Tab 4 - PCA
      pca <- prcomp(t(filtered_genes()[-1]), scale. = TRUE)
        # Store PCA results in a reactive value
      pca_results <- reactiveValues(
        pca_object = pca,
        explained_variance = round(100 * pca$sdev^2 / sum(pca$sdev^2), 2) # % variance explained
        )
        
        # Populate selectInput choices with available principal components
        updateSelectInput(
          session,
          inputId = "pc_x",
          choices = paste0("PC", 1:ncol(pca$x)),
          selected = "PC1"
        )
        updateSelectInput(
          session,
          inputId = "pc_y",
          choices = paste0("PC", 1:ncol(pca$x)),
          selected = "PC2"
        )
      
      # Scatter Plot of Selected Principal Components
        output$pca_scatter <- renderPlot({
          req(input$pc_x, input$pc_y, filtered_genes())
          
          pca <- pca_results$pca_object
          explained <- pca_results$explained_variance
          
          pc_x <- as.numeric(gsub("PC", "", input$pc_x))
          pc_y <- as.numeric(gsub("PC", "", input$pc_y))
          
          df <- data.frame(
            X = pca$x[, pc_x],
            Y = pca$x[, pc_y],
            Sample = rownames(pca$x)
          )
          
          ggplot(df, aes(x = X, y = Y, label = Sample)) +
            geom_point(color = "blue", size = 3) +
            geom_text_repel(aes(label = Sample)) +  # Use geom_text_repel
            labs(
              title = paste("PCA Scatter Plot:", input$pc_x, "vs", input$pc_y),
              x = paste0(input$pc_x, " (", explained[pc_x], "% Variance)"),
              y = paste0(input$pc_y, " (", explained[pc_y], "% Variance)")
            ) +
            theme_minimal()
        })
        
      
      # Beeswarm Plot of Top N Principal Components
      observeEvent(input$plot_beeswarm, {
        output$pca_beeswarm <- renderPlot({
          req(filtered_genes(), input$top_n)
          
          pca <- pca_results$pca_object
          explained <- pca_results$explained_variance
          top_n <- min(input$top_n, ncol(pca$x))  # Ensure top_n doesn't exceed available PCs
          
          df <- data.frame(
            PC = rep(paste0("PC", 1:top_n), each = nrow(pca$x)),
            Value = as.vector(pca$x[, 1:top_n])
          )
          
          ggplot(df, aes(x = PC, y = Value)) +
            geom_boxplot(outlier.color = "red", fill = "lightblue", alpha = 0.7) +
            labs(
              title = "Beeswarm Plot of Top Principal Components",
              x = "Principal Component",
              y = "PC Value"
            ) +
            theme_minimal()
        })
      })
      
      
    }, error = function(e) {
      output$error_message_series2 <- renderText("Oops! Wrong file. Please check the file.")
    })
    })
  
  
  
  ## Third Tab Output - differential expression
  deseq_results <- reactiveVal(NULL)
  observeEvent(input$run_analysis, {
    # Run DESeq and display results
    output$deseq_results <- renderDT({
      req(filtered_genes(), metadata_table(), input$design_factor)
      
      deseq_res <- run_deseq(filtered_genes(), metadata_table(), input$design_factor)
      
      deseq_results(deseq_res)  # Update reactiveVal
      
      datatable(deseq_results(),
                options = list(
                  pageLength = 8,
                  lengthMenu = c(5, 10, 20),
                  autoWidth = TRUE,
                  dom = 'Bfrtip',
                  buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
                ),
                rownames = TRUE,
                filter = 'top')
    })
    
    # Create labeled results reactively
    labeled_results <- reactive({
      req(deseq_results(), input$padj_threshold)  # Ensure inputs are valid
      label_res(deseq_results(), input$padj_threshold)
    })
    
    # Plot DESeq Results
    output$pvalue_hist <- renderPlot({
      req(labeled_results())  # Ensure labeled_results is ready
      plot_pvals(labeled_results())
    })
    
    output$log2fc_hist <- renderPlot({  # Fixed typo
      req(labeled_results())
      plot_log2fc(labeled_results(), input$padj_threshold)
    })
    
    output$volcano_scatter <- renderPlot({
      req(labeled_results())
      plot_volcano(labeled_results())
    })
  })
  
  
  ## Fourth Tab Output - Network
  data <- reactiveVal(NULL)
  gene_cor <- reactiveVal(NULL)
  graph_obj <- reactiveVal(NULL)
  # Load data
  observeEvent(input$datafile, {
    req(input$datafile)
    data(read.csv(input$datafile$datapath, row.names = 1))
  })
  
  # Subset genes and compute pairwise correlations
  observeEvent(input$gene_list, {
    req(data())
    gene_names <- unlist(strsplit(input$gene_list, "\n"))
    gene_names <- trimws(gene_names)
    
    # Check if all genes are present in the data
    missing_genes <- setdiff(gene_names, rownames(data()))
    if (length(missing_genes) > 0) {
      showNotification(
        paste("Missing genes:", paste(missing_genes, collapse = ", ")),
        type = "error"
      )
    }
    
    # Subset and compute correlations
    valid_genes <- intersect(gene_names, rownames(data()))
    subset_data <- data()[valid_genes, ]
    correlations <- cor(t(subset_data), method = "pearson")
    gene_cor(correlations)
    
    # Update vertex dropdown choices
    updateSelectizeInput(session, "vertex1", choices = valid_genes, server = TRUE)
    updateSelectizeInput(session, "vertex2", choices = valid_genes, server = TRUE)
  })
  
  # Create heatmap
  output$heatmap <- renderPlot({
    req(gene_cor())
    pheatmap(
      gene_cor(),
      main = "Clustered Heatmap of Gene Correlations",
      xlab = "Genes",
      ylab = "Genes",
      cluster_rows = TRUE,
      cluster_cols = TRUE
    )
  })
  
  # Create Correlation network
  observeEvent(input$cor_threshold, {
    req(gene_cor())  # Ensure gene correlation matrix is available
    threshold <- input$cor_threshold
    
    # Create graph edges based on the threshold
    edges <- which(gene_cor() > threshold, arr.ind = TRUE)
    edges <- edges[edges[, 1] != edges[, 2], ]  # Remove self-loops
    
    # Check if there are any edges above the threshold
    if (nrow(edges) == 0) {
      # If no edges exist, create a graph with isolated vertices
      graph <- make_empty_graph(n = nrow(gene_cor()), directed = FALSE)
      V(graph)$name <- rownames(gene_cor())  # Assign names to vertices
    } else {
      # If edges exist, create graph from the edge list
      graph <- graph_from_edgelist(as.matrix(edges), directed = FALSE)
      
      # Assign edge weights
      E(graph)$weight <- gene_cor()[edges]
      
      # Ensure all genes are represented as vertices (even isolated ones)
      vertex_names <- rownames(gene_cor())
      if (length(vertex_names) != vcount(graph)) {
        # Add missing vertices if any are isolated
        missing_vertices <- setdiff(vertex_names, V(graph)$name)
        graph <- add_vertices(graph, length(missing_vertices), name = missing_vertices)
      }
      
      # Assign vertex names
      V(graph)$name <- vertex_names
    }
    
    # Update graph object
    graph_obj(graph)
  })
  
  output$network_plot <- renderPlot({
    req(graph_obj())
    graph <- graph_obj()
    
    plot(
      graph,
      vertex.label = V(graph)$name,
      vertex.label.cex = 0.8,           # Adjust label size
      vertex.label.color = "black",    # Set label color
      vertex.label.dist = 2,           # Push labels away from nodes
      vertex.color = rainbow(vcount(graph)),
      vertex.size = 15,
      edge.width = E(graph)$weight * 5,
      layout = layout_with_fr(graph)   # Fruchterman-Reingold layout for better spacing
    )
  })
  
  # Compute network metrics
  output$network_metrics <- renderDataTable({
    req(graph_obj())
    graph <- graph_obj()
    metrics <- data.frame(
      Degree = degree(graph),
      Closeness = closeness(graph),
      Betweenness = betweenness(graph)
    ) %>%
      rownames_to_column(var = 'Gene')
    
    datatable(metrics, options = list(pageLength = 10))
  })
  
  # Shortest path
  observeEvent(input$shortest_path_btn, {
    req(graph_obj(), input$vertex1, input$vertex2)
    graph <- graph_obj()
    path <- shortest_paths(graph, from = input$vertex1, to = input$vertex2, output = "vpath")
    
    output$shortest_path <- renderText({
      if (length(path$vpath[[1]]) > 0) {
        paste("Shortest path:", paste(V(graph)$name[path$vpath[[1]]], collapse = " -> "))
      } else {
        "No path exists between the selected vertices."
      }
    })
  })
}


# Run the application 
shinyApp(ui = ui, server = server)
