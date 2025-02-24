library(ggplot2)
library(plotly)
library(shiny)
library(dplyr)
library(DESeq2)

options(shiny.maxRequestSize = 100 * 1024^2)  # Increase upload limit to 100MB

# Define UI
ui <- fluidPage(
  titlePanel("good.genes"),
  sidebarLayout(
    sidebarPanel(
      fileInput("deseq_file", "Input CSV", accept = ".csv"),
      sliderInput("x_intercept", "Fold Change:", min = 0, max = 50, value = 0, step = 1),
      sliderInput("parabola_b", "Significance:", min = 0, max = 350, value = 0, step = 1),
      sliderInput("logBaseMean", "Log10(BaseMean):", min = 0, max = 5.5, value = c(0, 5.5), step = 0.1),
      ),
  mainPanel(
      plotlyOutput("volcanoPlot", height = "700px")
    )
  )
)

# Define Server
server <- function(input, output, session) {
  observe({
    req(input$parabola_a > 0, input$parabola_b > 0)
    x_left <- -sqrt(input$parabola_b / input$parabola_a)
    x_right <- sqrt(input$parabola_b / input$parabola_a)
    updateSliderInput(session, "parabola_a", min = x_left, max = x_right)
  })
  
  data_reactive <- reactive({
    req(input$deseq_file)
    
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
    
    # Handle padj values of zero by setting them to 1e-310 to avoid zero rounding error
    deseq_data$padj[deseq_data$padj == 0] <- 1e-310
    
    # Transform padj to -log10 scale, ensuring no negative values
    deseq_data <- deseq_data %>% mutate(
      negLog10P = ifelse(!is.na(padj) & padj > 0, -log10(padj), NA),
      logBaseMean = log10(baseMean + 1)  # Avoids log(0)
    )
    
    return(deseq_data)
  })
  
  filtered_data <- reactive({
  req(data_reactive())
  df <- data_reactive()

  df_filtered <- df %>%
    filter(
      logBaseMean >= input$logBaseMean[1],
      logBaseMean <= input$logBaseMean[2]
    ) %>%
    filter(
      negLog10P >= (-(input$parabola_b / (input$x_intercept^2)) * (log2FoldChange^2) + input$parabola_b)
    )

  return(df_filtered)
  })
  
  output$volcanoPlot <- renderPlotly({
    req(filtered_data())
    gg <- ggplot(filtered_data(), aes(x = log2FoldChange, y = negLog10P, size = baseMean, text = paste0("Gene: ", gene))) +
      geom_point(alpha = 0.3, color = "blue") +
      scale_size_continuous(range = c(3, 15), limits = range(data_reactive()$baseMean, na.rm = TRUE)) +
      labs(x = "Log2 Fold Change", y = "-Log10(Adjusted P-value)", title = "Interactive Volcano Bubble Plot") +
      theme_minimal()
    
    ggplotly(gg, tooltip = "text") %>%
      layout(yaxis = list(range = c(min(data_reactive()$negLog10P, na.rm = TRUE) - 5, max(data_reactive()$negLog10P, na.rm = TRUE) + 5)), 
             xaxis = list(range = c(min(data_reactive()$log2FoldChange, na.rm = TRUE) - 5, max(data_reactive()$log2FoldChange, na.rm = TRUE) + 5)))
  })
}

# Run App
shinyApp(ui = ui, server = server)
