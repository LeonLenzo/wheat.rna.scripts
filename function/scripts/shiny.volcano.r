library(ggplot2)
library(plotly)
library(shiny)
library(dplyr)
library(DESeq2)

options(shiny.maxRequestSize = 100 * 1024^2)  # Increase upload limit to 100MB

# Define UI
ui <- fluidPage(
  titlePanel("Interactive Volcano Bubble Plot"),
  sidebarLayout(
    sidebarPanel(
      fileInput("deseq_file", "Upload DESeq2 Results (CSV)", accept = ".csv"),
      fileInput("go_file", "Upload GO Annotations (CSV)", accept = ".csv"),
      sliderInput("log2fc", "Log2 Fold Change Threshold:", min = 0, max = 50, value = 2),
      sliderInput("pval", "Adjusted P-value Threshold (-log10):", min = 0, max = 310, value = 0),
      sliderInput("basemean", "BaseMean Range:", min = 0, max = 10000, value = c(0, 10000)),
      selectizeInput("go_category", "GO Category:", choices = NULL, multiple = TRUE),
      selectizeInput("go_term", "GO Term Description:", choices = NULL, multiple = TRUE, options = list(maxOptions = 500))
    ),
    mainPanel(
      plotlyOutput("volcanoPlot", height = "700px")
    )
  )
)

# Define Server
server <- function(input, output, session) {
  
  data_reactive <- reactive({
    req(input$deseq_file, input$go_file)
    
    # Read DESeq2 results
    deseq_data <- read.csv(input$deseq_file$datapath, stringsAsFactors = FALSE)
    
    if (!"gene" %in% colnames(deseq_data)) {
      stop("DESeq2 results file must contain a 'gene' column.")
    }
    
    # Ensure gene column is character type
    deseq_data$gene <- as.character(deseq_data$gene)
    
    if (nrow(deseq_data) == 0) {
      return(data.frame())
    }
    
    # Handle padj values of zero by setting them to 1e-310
    deseq_data$padj[deseq_data$padj == 0] <- 1e-310
    
    # Transform padj to -log10 scale, ensuring no negative values
    deseq_data <- deseq_data %>% mutate(negLog10P = ifelse(!is.na(padj) & padj > 0, -log10(padj), NA))
    
    # Read GO annotations and rename columns
    go_data <- read.csv(input$go_file$datapath, stringsAsFactors = FALSE)
    colnames(go_data) <- c("gene", "GO_term", "GO_description", "GO_category")
    go_data$gene <- as.character(go_data$gene)  # Ensure gene column is character
    
    if (nrow(go_data) == 0) {
      return(data.frame())
    }
    
    if (!all(c("gene", "GO_category", "GO_term", "GO_description") %in% colnames(go_data))) {
      stop("Uploaded GO annotations file must contain Gene_ID, GO_Terms, GO_Description, and GO_Category columns.")
    }
    
    # Merge GO and DESeq2 data
    merged_data <- left_join(deseq_data, go_data, by = "gene")
    
    if (nrow(merged_data) > 0) {
      updateSelectizeInput(session, "go_category", choices = unique(na.omit(merged_data$GO_category)), server = TRUE)
      updateSelectizeInput(session, "go_term", choices = unique(na.omit(merged_data$GO_description)), server = TRUE)
    }
    
    return(merged_data)
  })
  
  filtered_data <- reactive({
    req(data_reactive())
    df <- data_reactive()
    
    go_category_filter <- if (length(input$go_category) == 0) unique(df$GO_category) else input$go_category
    go_term_filter <- if (length(input$go_term) == 0) unique(df$GO_description) else input$go_term

    df_filtered <- df %>%
  group_by(gene, log2FoldChange, negLog10P, baseMean) %>%
  summarize(
    GO_category = ifelse(all(is.na(GO_category)), "None", paste(unique(na.omit(GO_category)), collapse = ", ")),
    GO_description = ifelse(all(is.na(GO_description)), "None", paste(unique(na.omit(GO_description)), collapse = "; ")),
    .groups = "drop"
  ) %>%
  group_by(gene, log2FoldChange, negLog10P, baseMean) %>%
  summarize(
    GO_category = paste(unique(GO_category), collapse = ", "),
    GO_description = paste(unique(GO_description), collapse = "; "),
    .groups = "drop"
  ) %>%
      filter(
  (GO_category %in% go_category_filter | "None" %in% go_category_filter),
  (GO_description %in% go_term_filter | "None" %in% go_term_filter),
        abs(log2FoldChange) >= input$log2fc | input$log2fc == 0,
        !is.na(negLog10P) & negLog10P >= input$pval,
        baseMean >= input$basemean[1] & baseMean <= input$basemean[2],
        GO_category %in% go_category_filter,
        GO_description %in% go_term_filter
      )
    
    return(df_filtered)
  })
  
  output$volcanoPlot <- renderPlotly({
    req(filtered_data())
    gg <- ggplot(filtered_data(), aes(x = log2FoldChange, y = negLog10P, size = baseMean, color = GO_category, text = paste0("Gene: ", gene, 
              "<br>GO Categories: ", GO_category, 
              "<br>GO Terms: ", GO_description))) +
      geom_point(alpha = 0.7) +
      scale_size(range = c(3, 15)) +
      labs(x = "Log2 Fold Change", y = "-Log10(Adjusted P-value)", title = "Interactive Volcano Bubble Plot") +
      theme_minimal()
    
    ggplotly(gg, tooltip = "text")
  })
}

# Run App
shinyApp(ui = ui, server = server)
