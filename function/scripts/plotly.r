# Load required libraries
library(ggplot2)
library(dplyr)
library(readr)
library(plotly)
library(htmlwidgets)

# Function to generate and save an interactive Plotly bubble chart
generate_interactive_bubble_chart <- function(csv_file, output_html) {
  # Read the CSV file
  go_results <- read_csv(csv_file, show_col_types = FALSE)

  # Ensure required columns exist
  required_cols <- c("GO_Description", "p.adjust", "Count", "Sum_log2FoldChange", "GO_Category")
  if (!all(required_cols %in% colnames(go_results))) {
    message("Skipping file due to missing columns: ", csv_file)
    return(NULL)
  }

  # Convert necessary columns to numeric and clean data
  go_results <- go_results %>%
    mutate(
      Sum_log2FoldChange = as.numeric(Sum_log2FoldChange),
      p.adjust = as.numeric(p.adjust),
      GO_Category = as.factor(GO_Category),
      Tooltip = paste("GO Term:", GO_Description, 
                      "<br>Category:", GO_Category,
                      "<br>-log10 p-value:", round(-log10(p.adjust), 3),
                      "<br>Fold Change:", round(Sum_log2FoldChange, 3))
    ) %>%
    filter(!is.na(Sum_log2FoldChange) & !is.na(p.adjust))

  # Create a ggplot bubble chart
  plot <- ggplot(go_results, aes(x = Sum_log2FoldChange, y = -log10(p.adjust), text = Tooltip)) +
    geom_point(aes(size = Count, fill = GO_Category), 
               shape = 21, color = "transparent", alpha = 0.8) +
    scale_fill_manual(values = c("Biological Process" = "#30343F", 
                                 "Molecular Function" = "#656839", 
                                 "Cellular Component" = "#89023E")) +
    scale_size(range = c(3, 10)) +
    theme_minimal() +
    labs(
      title = paste("GO Enrichment:", tools::file_path_sans_ext(basename(csv_file))),
      x = "Sum of log2 Fold Change",
      y = "-log10 Adjusted p-value",
      fill = "GO Category",
      size = "Gene Count"
    )

  # Convert ggplot to interactive Plotly chart
  interactive_plot <- ggplotly(plot, tooltip = "text")

  # Save as an interactive HTML file
  saveWidget(interactive_plot, output_html, selfcontained = TRUE)

  message("✅ Saved:", output_html)
}

# Function to process all CSV files in a directory
process_directory <- function(input_dir) {
  # Get all CSV files in the directory
  csv_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

  if (length(csv_files) == 0) {
    stop("No CSV files found in the directory:", input_dir)
  }

  # Process each CSV file and generate an HTML output
  for (csv_file in csv_files) {
    output_html <- file.path(input_dir, paste0(tools::file_path_sans_ext(basename(csv_file)), "_BubbleChart.html"))
    tryCatch({
      generate_interactive_bubble_chart(csv_file, output_html)
    }, error = function(e) {
      message("❌ Error processing file: ", csv_file, "\n", e$message)
    })
  }

  message("🎉 Processing complete! HTML files saved in:", input_dir)
}

# Example Usage:
# Replace `"your_directory_path"` with the actual folder containing CSV files
input_dir <- "GO"  # Change this to your actual directory
process_directory(input_dir)
