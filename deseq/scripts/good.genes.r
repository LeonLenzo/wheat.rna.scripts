library(shiny)
library(plotly)
library(dplyr)
library(ggplot2)
library(DT)
library(shinyBS)
library(rsconnect)


ui <- fluidPage(
  titlePanel("good.genes"),
  sidebarLayout(
    sidebarPanel(
      h3("Controls"),
      textInput("search_input", "Search Gene or Keyword:", placeholder = "Type a gene name or keyword..."),
      checkboxInput("show_unselected", "Show Unselected Genes", value = TRUE),  # <-- New checkbox
      bsCollapse(
        id = "sidebar_panels",
        open = "Inputs",
        bsCollapsePanel(
          "Inputs",
          fileInput("deseq_files", "Upload Expression Data (CSV)", accept = ".csv", multiple = TRUE),
          selectInput("dataset_selector", "Select Dataset:", choices = NULL),
          fileInput("annotation_file", "Upload Annotation Data (CSV)", accept = ".csv"),
          style = "primary"
        ),
        bsCollapsePanel(
          "Filters",
          numericInput("x_intercept", "Fold Change Min:", value = 1, min = 0, max = 50, step = 1),
          numericInput("parabola_b", "Significance Min:", value = 1, min = 0, max = 350, step = 1),
          numericInput("logBaseMean_min", "Minimum Expression:", value = 0, min = 0, max = 5.5, step = 0.1),
          style = "success"
        ),
        bsCollapsePanel(
          "Selected Genes",
          DTOutput("selected_genes_table"),
          downloadButton("download_genes", "Download Table"),
          style = "info"
        )
      )
    ),
    mainPanel(
      plotlyOutput("volcanoPlot", height = "700px"),
      hr(),
      h4("Selected Gene:"),
      uiOutput("gene_info")
    )
  )
)


server <- function(input, output, session) {
  # Store multiple datasets
  datasets <- reactiveValues(data = list())

  # Load and process uploaded DESeq2 datasets
  observeEvent(input$deseq_files, {
    req(input$deseq_files)
    uploaded_files <- input$deseq_files$datapath
    file_names <- input$deseq_files$name

    datasets$data <- lapply(uploaded_files, function(path) {
      df <- read.csv(path, stringsAsFactors = FALSE)
      if (!"gene" %in% colnames(df)) {
        stop("DESeq2 results file must contain a 'gene' column.")
      }
      df$gene <- as.character(df$gene)
      df$padj[df$padj == 0] <- 1e-310  # Avoid log(0) issues
      df <- df %>% mutate(
        negLog10P = ifelse(!is.na(padj) & padj > 0, -log10(padj), NA),
        logBaseMean = log10(baseMean + 1)
      )
      return(df)
    })
    names(datasets$data) <- file_names  # Name datasets by their file names

    # Update dataset selection dropdown
    updateSelectInput(session, "dataset_selector", choices = file_names, selected = file_names[1])
  })

  # Reactive expression to load the annotation file
  annotation_data <- reactive({
    req(input$annotation_file)
    read.csv(input$annotation_file$datapath, stringsAsFactors = FALSE)
  })

  # Cached joined data (using annotation only once)
  joined_data <- reactive({
    req(input$dataset_selector, datasets$data)
    df <- datasets$data[[input$dataset_selector]]
    if (!is.null(input$annotation_file)) {
      annot <- annotation_data()
      df <- dplyr::left_join(df, annot, by = c("gene" = "prot_ID"))
    }
    df
  })

  # Reactive storage for selected gene and GO term
  selected <- reactiveValues(gene = NULL, go_term = NULL)

  filtered_data <- reactive({
    req(joined_data())
    df <- joined_data() %>%
      filter(logBaseMean >= input$logBaseMean_min) %>%
      filter(negLog10P >= (-(input$parabola_b / (input$x_intercept^2)) * (log2FoldChange^2) + input$parabola_b))
    df
  })

  observeEvent(event_data("plotly_click"), {
    click_data <- event_data("plotly_click")
    if (!is.null(click_data) && "key" %in% names(click_data)) {
      selected$gene <- click_data$key
      # Do not clear the GO term so that the set remains selected.
      # selected$go_term <- NULL
    }
  })


  observeEvent(event_data("plotly_doubleclick"), {
    selected$gene <- NULL
    selected$go_term <- NULL
  })

  # Track the clicked GO term
  observeEvent(input$go_term_click, {
    if (!is.null(input$go_term_click)) {
      selected$go_term <- input$go_term_click
    }
  })

  # Find genes that share the selected GO term
  shared_genes <- reactive({
    req(filtered_data(), selected$go_term)

    shared <- filtered_data() %>%
      filter(
        grepl(selected$go_term, Molecular.Function, fixed = TRUE) |
        grepl(selected$go_term, Biological.Process, fixed = TRUE) |
        grepl(selected$go_term, Cellular.Component, fixed = TRUE)
      )
    return(shared)
  })


output$volcanoPlot <- renderPlotly({
  req(filtered_data())
  
  plot_data <- filtered_data()
  
  # Compute fixed axis ranges based on all filtered data
  x_range <- range(plot_data$log2FoldChange, na.rm = TRUE)
  y_range <- range(plot_data$negLog10P, na.rm = TRUE)
  
  # Compute margins (5% of the range, adjust as needed)
  margin_x <- diff(x_range) * 0.05
  margin_y <- diff(y_range) * 0.05
  
  # Compute global bubble size limits from the complete data
  global_minBM <- min(plot_data$baseMean, na.rm = TRUE)
  global_maxBM <- max(plot_data$baseMean, na.rm = TRUE)
  
  highlight_gene <- selected$gene
  highlight_go_term <- selected$go_term
  search_keyword <- input$search_input
  
  # Default to no GO term highlight
  go_highlight_color <- NULL
  go_shared_points <- NULL
  if (!is.null(highlight_go_term)) {
    go_shared_points <- shared_genes()
    if (nrow(go_shared_points) > 0) {
      if (any(grepl(highlight_go_term, go_shared_points$Molecular.Function, fixed = TRUE))) {
        go_highlight_color <- "#9b59b6"  # Purple for Molecular Function
      } else if (any(grepl(highlight_go_term, go_shared_points$Biological.Process, fixed = TRUE))) {
        go_highlight_color <- "#2ecc71"  # Green for Biological Process
      } else if (any(grepl(highlight_go_term, go_shared_points$Cellular.Component, fixed = TRUE))) {
        go_highlight_color <- "#3498db"  # Blue for Cellular Component
      }
    }
  }
  
  # Identify the selected gene
  selected_point <- NULL
  if (!is.null(highlight_gene) && highlight_gene %in% plot_data$gene) {
    selected_point <- plot_data %>% filter(gene == highlight_gene)
  }
  
  # Identify genes matching the search
  search_matches <- NULL
  if (!is.null(search_keyword) && search_keyword != "") {
    search_matches <- plot_data %>% filter(
      grepl(search_keyword, gene, ignore.case = TRUE) | 
      grepl(search_keyword, Description, ignore.case = TRUE)
    )
  }
  
  # Build the ggplot object with conditional layers.
  # Notice we now add explicit limits to scale_size_continuous.
  gg <- ggplot() +
    { if (input$show_unselected) 
        geom_point(data = plot_data, aes(
          x = log2FoldChange, 
          y = negLog10P,
          size = baseMean,
          key = gene
        ), alpha = 0.3, color = "#D9D9D9")
    } +
    { if (!is.null(go_shared_points) && nrow(go_shared_points) > 0 && !is.null(go_highlight_color)) 
        geom_point(data = go_shared_points, aes(
          x = log2FoldChange, 
          y = negLog10P,
          size = baseMean,
          key = gene
        ), alpha = 0.6, color = go_highlight_color)
    } +
    { if (!is.null(search_matches) && nrow(search_matches) > 0 && is.null(go_highlight_color)) 
        geom_point(data = search_matches, aes(
          x = log2FoldChange, 
          y = negLog10P,
          size = baseMean,
          key = gene
        ), alpha = 0.7, color = "#FF8000")
    } +
    { if (!is.null(selected_point) && nrow(selected_point) > 0) 
        geom_point(data = selected_point, aes(
          x = log2FoldChange, 
          y = negLog10P,
          size = baseMean,
          key = gene
        ), alpha = 0.9, color = "#A50021")
    } +
    scale_size_continuous(range = c(3, 15), limits = c(global_minBM, global_maxBM)) +
    labs(x = "Log2 Fold Change", y = "-Log10(Adjusted P-value)") +
    theme_minimal()
  
  # Convert the ggplot object to a Plotly object
  p <- ggplotly(gg, tooltip = "none") %>% 
    event_register("plotly_click") %>% 
    layout(
      hovermode = 'closest',
      xaxis = list(title = "Log2 Fold Change", automargin = TRUE,
                   range = c(x_range[1] - margin_x, x_range[2] + margin_x)),
      yaxis = list(title = "-Log10(Adjusted P-value)", automargin = TRUE,
                   range = c(y_range[1] - margin_y, y_range[2] + margin_y))
    )
  
  # Convert to WebGL and remove unsupported attributes
  p_webgl <- plotly::toWebGL(p)
  p_webgl$x$data <- lapply(p_webgl$x$data, function(trace) {
    if ("hoveron" %in% names(trace)) trace$hoveron <- NULL
    trace
  })
  
  p_webgl
})



  # Render gene information with clickable GO terms
  output$gene_info <- renderUI({
    req(selected$gene)
    
    gene_data <- filtered_data() %>% filter(gene == selected$gene)
    if (nrow(gene_data) == 0) return(HTML("<i>No data for selected gene</i>"))
    
    # Helper function to create clickable GO term links
    create_go_link <- function(go_terms) {
      if (is.null(go_terms) || is.na(go_terms) || go_terms == "") {
        return("N/A")
      }
      
      # Split GO terms by common delimiters (`,`, `;`, or `|`)
      go_terms_list <- unlist(strsplit(go_terms, "[;]"))
      go_terms_list <- trimws(go_terms_list)
      
      # Create individual clickable links for each GO term
      links <- sapply(go_terms_list, function(term) {
        sprintf('<a href="#" onclick="Shiny.setInputValue(\'go_term_click\', \'%s\', {priority: \'event\'}); return false;">%s</a>', term, term)
      })
      return(paste(links, collapse = ", "))
    }
    
    HTML(paste(
      "<b>Gene:</b> ", gene_data$gene, "<br>",
      "<b>Log2FC:</b> ", round(gene_data$log2FoldChange, 2), "<br>",
      "<b>-Log10(padj):</b> ", round(gene_data$negLog10P, 2), "<br>",
      "<b>Base Mean:</b> ", round(gene_data$baseMean, 2), "<br>",
      "<b>Description:</b> ", ifelse(is.na(gene_data$Description), "N/A", gene_data$Description), "<br>",
      "<span style='color:purple;'><b>Molecular Function:</b></span> ", create_go_link(gene_data$Molecular.Function), "<br>",
      "<span style='color:green;'><b>Biological Process:</b></span> ", create_go_link(gene_data$Biological.Process), "<br>",
      "<span style='color:blue;'><b>Cellular Component:</b></span> ", create_go_link(gene_data$Cellular.Component)
    ))
  })

  output$selected_genes_table <- renderDT({
    gene_list <- filtered_data()

    # If a GO term is selected, filter by it first
    if (!is.null(selected$go_term)) {
      gene_list <- shared_genes()
    }

    # Apply master search filter (gene name OR description keyword)
    if (!is.null(input$search_input) && input$search_input != "") {
      gene_list <- gene_list %>% 
        filter(
          grepl(input$search_input, gene, ignore.case = TRUE) | 
          (if ("Description" %in% colnames(gene_list)) grepl(input$search_input, Description, ignore.case = TRUE) else FALSE)
        )
    }

    if (nrow(gene_list) == 0) return(NULL)

    # Make gene names clickable & round values
    gene_list <- gene_list %>% 
      mutate(
        Gene = sprintf('<a href="#" onclick="Shiny.setInputValue(\'go_gene_click\', \'%s\', {priority: \'event\'}); return false;">%s</a>', gene, gene),
        log2FoldChange = round(log2FoldChange, 0),
        negLog10P = round(negLog10P, 0),
        baseMean = round(baseMean, 0)
      )

    # Only include Description if it exists
    table_data <- if ("Description" %in% colnames(gene_list)) {
      select(gene_list, Gene, Description, log2FoldChange, negLog10P, baseMean)
    } else {
      select(gene_list, Gene, log2FoldChange, negLog10P, baseMean)
    }

    colnames(table_data) <- c("Gene", "Description", "Log2 Fold Change", "-Log10(p-adj)", "Base Mean")[1:ncol(table_data)]

    DT::datatable(
      table_data,
      escape = FALSE,
      options = list(pageLength = 10, dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel'))
    )
  }, server = TRUE)


  # Track when a gene in the table is clicked
  observeEvent(input$go_gene_click, {
    selected$gene <- input$go_gene_click
    #selected$go_term <- NULL  # Reset GO term when selecting a gene
  })

  # Download selected genes table as CSV
  output$download_genes <- downloadHandler(
    filename = function() {
      paste0("selected_genes_", Sys.Date(), ".csv")
    },
    content = function(file) {
      if (!is.null(selected$go_term)) {
        go_shared <- shared_genes()
        if (nrow(go_shared) > 0) {
          write.csv(go_shared %>% select(gene, Description, log2FoldChange, negLog10P, baseMean), file, row.names = FALSE)
        }
      } else if (!is.null(selected$gene)) {
        gene_data <- filtered_data() %>% filter(gene == selected$gene)
        if (nrow(gene_data) > 0) {
          write.csv(gene_data %>% select(gene, Description, log2FoldChange, negLog10P, baseMean), file, row.names = FALSE)
        }
      }
    }  
  )
} 

shinyApp(ui = ui, server = server)