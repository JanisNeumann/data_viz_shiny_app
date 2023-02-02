#### Dependencies ####

library(shiny)
library(ggplot2)
library(dplyr)
library(disgenet2r)
library(rbioapi)



#### UI ####
ui <- fluidPage(
  
  # Application title
  titlePanel("Data Visualization using DisGeNet & Reactome"),
  sidebarLayout(
    sidebarPanel(
      # Define inputs as two .csv files
      fileInput("gene_data", "Choose CSV file containing data",
                accept = c(".csv")
      ),
      textInput("dgn_acct", "DisGeNet: Email", value = ""),
      textInput("dgn_pwd", "DisGeNet: Password", value = "")
    ),
    
    # Define UMAP plot as output
    mainPanel(
      tabsetPanel(
        tabPanel("Overview", tableOutput("count_table")),
        tabPanel("Volcano Plot", plotOutput("volcano_plot", height = 1000, width = 1100)),
        tabPanel("Retained Identifiers", tableOutput("retained")),
        tabPanel("Gene-Disease Network", plotOutput("dgn_net", height = 1000, width = 1100)),
        tabPanel("Gene-Disease Heat Map", plotOutput("dgn_hm", height = 1000, width = 1100)),
        tabPanel("Enriched Pathways Table", tableOutput("enriched_table")),
        tabPanel("Enriched Pathways Plot", plotOutput("enriched_plot", height = 1000, width = 750))
      )
    )
  )
)



#### Server ####
server <- function(input, output, session) {
  
  #### Constants ####
  
  padjust_threshold <- 0.05
  FC_threshold <- 1
  reactome_padjust_threshold <- padjust_threshold
  
  
  
  #### Input-dependent variables ####
  
  # Reactive expression that reads the file
  gene_data <- reactive({
    req(input$gene_data)
    
    read.csv(input$gene_data$datapath, header = TRUE, sep = ",")
    
  })
  
  # Reformat some data for the volcano plot and data overview.
  volcano_data <- reactive({
    vpd <- gene_data()
    if (is.null(vpd)) {
      return(NULL)
    }
    
    # Record down- and upregulation according to p-value & log fold change
    tmp_reg <- rep("Abs. log FC < 1", nrow(vpd))
    tmp_reg[vpd$adj.pValue <= padjust_threshold & vpd$logFC <= -FC_threshold] <- "Downregulated"
    tmp_reg[vpd$adj.pValue <= padjust_threshold & vpd$logFC >= FC_threshold] <- "Upregulated"
    tmp_reg <- factor(tmp_reg, levels = c("Abs. log FC < 1", "Downregulated", "Upregulated"))
    
    vpd$Regulation <- tmp_reg
    
    vpd$Significant <- as.factor(ifelse(vpd$adj.pValue <= padjust_threshold, "Yes", "No"))
    
    vpd
    
  })
  
  # Subset the data to retain genes/proteins with adj. p-value <= 0.05 & abs. log FC >= 1
  retained_data <- reactive({
    rd <- gene_data()
    if (is.null(rd)) {
      return(NULL)
    }
    
    rd[rd$adj.pValue <= padjust_threshold & abs(rd$logFC) >= FC_threshold, ]
    
  })
  
  # DisGeNet Setup with user-provided account information: Get the API key.
  dgn_api <- reactive({
    req(retained_data,
        input$dgn_acct,
        input$dgn_pwd
    )
    
    dgn_acct <- input$dgn_acct
    dgn_pwd <- input$dgn_pwd
    
    # Retrieve API key.
    dgn_api <- get_disgenet_api_key(
      email = dgn_acct, 
      password = dgn_pwd)
    
    dgn_api
    
  })
  
  # DisGeNet Setup with user-provided account information: Use API key to get information about gene2disease data.
  disgenet_info <- reactive({
    req(dgn_api)
    
    disgenet2r::gene2disease(retained_data()$Symbol,
                             vocabulary = "HGNC",
                             database = "CURATED",
                             score = c(0.6, 1),
                             api_key = dgn_api(),
                             verbose = FALSE,
                             warnings = FALSE
    )
  })
  
  # Use enrichr integration of rbioapi package to get Reactome enrichment data.
  enrichr_data <- reactive({
    req(retained_data)
    
    # Get enrichment terms and remove weird columns relating to older p-value calculations.
    # Also filter out insignificant terms.
    rba_enrichr(gene_list = retained_data()$Symbol,
                gene_set_library = "Reactome_2022") %>%
      select_if(!names(.) %in% c("Old.P.value", "Old.Adjusted.P.value")) %>%
      filter(Adjusted.P.value <= reactome_padjust_threshold)
    
  })
  
  # Terms sorted by adj. p-val, Overlap as a number, and adj. p-val for enrichment plot.
  enriched_plot_data <- reactive({
    req(enrichr_data)
    
    enriched_plot_data <- data.frame(Term = factor(enrichr_data()$Term, levels = enrichr_data()$Term[order(enrichr_data()$Adjusted.P.value, decreasing = TRUE)]),
                                     Overlap = as.numeric(gsub("/.*", "", enrichr_data()$Overlap))/as.numeric(gsub("^.*/", "", enrichr_data()$Overlap)),
                                     Adjusted_P_value = enrichr_data()$Adjusted.P.value)
    
  })
  
  
  
  #### Outputs ####
  
  # Tabular overview of significance, up- & downregulation etc.
  output$count_table <- renderTable({
    count_df <- volcano_data()
    if(is.null(count_df)){return()}
    
    data.frame("Observations" = c(nrow(count_df)),
               "Significant" = sum(count_df$Significant == "Yes"),
               "Downregulated" = sum(count_df$Regulation == "Downregulated"),
               "Upregulated" = sum(count_df$Regulation == "Upregulated"),
               row.names = NULL)
  })
  
  # Reactive expression for volcano plot.
  output$volcano_plot <- renderPlot({
    volcano_plot_df <- volcano_data()
    if(is.null(volcano_plot_df)){return()}
    
    volcano_plot <- ggplot(data = volcano_plot_df, aes(x = logFC, y = -log10(adj.pValue), col = Regulation, shape = Significant)) +
      geom_point() +
      theme_bw() +
      scale_shape_manual(values = c(1, 16)) +
      scale_color_manual(values = c("black", "firebrick3", "cyan3")) +
      geom_vline(xintercept = c(-FC_threshold, FC_threshold)) +
      geom_hline(yintercept = -log10(padjust_threshold))
    
    volcano_plot
  })
  
  # Tabular overview of all retained genes/proteins.
  output$retained <- renderTable({
    retained_df <- retained_data()
    if(is.null(retained_df)){return()}
    
    retained_df
    
  }, digits = 5)
  
  # Gene-disease association network according to DisGeNet.
  output$dgn_net <- renderPlot({
    
    dgn_d <- disgenet_info()
    
    disgenet2r::plot(dgn_d, class = "Network", prop = 20)
    
  })
  
  # Gene-disease association heat map according to DisGeNet.
  output$dgn_hm <- renderPlot({
    
    dgn_d <- disgenet_info()
    
    disgenet2r::plot(dgn_d, class = "DiseaseClass", prop = 3)
    
  })
  
  # Simply the table generated by Reactome enrichment analysis presented as a table.
  output$enriched_table <- renderTable({
    enriched_table <- enrichr_data()
    if(is.null(enriched_table)){
      return(NULL)
      }
    
    enriched_table
    
  }, digits = 5)
  
  # A ggplot scatter plot of Reactome enrichment analysis, with color corresponding to significance and point size to overlap.
  output$enriched_plot <- renderPlot({
    
    epd <- enriched_plot_data()
    
    ggplot(epd, aes(x = Adjusted_P_value, y = Term, color = Adjusted_P_value, size = Overlap)) +
      geom_point() +
      theme_bw() +
      xlab("Adjusted p-value")
    
  })
  
}



#### Run the app ####

shinyApp(ui, server)
