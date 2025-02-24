# Load required libraries
library(clusterProfiler)
library(ggplot2)
library(dplyr)
library(forcats)
library(readr)
library(ggrepel)  # For better label positioning
library(stringr)  # For text wrapping

# Function to generate a GO enrichment bubble chart
generate_bubble_chart <- function(csv_file, output_pdf) {
  # Read GO enrichment results from CSV
  go_results <- read_csv(csv_file, show_col_types = FALSE)
  
  # Check if required columns exist
  required_cols <- c("GO_Description", "p.adjust", "Count", "Sum_log2FoldChange", "GO_Category")
  if (!all(required_cols %in% colnames(go_results))) {
    stop(paste("Error: Missing required columns in", csv_file))
  }
  
  # Convert necessary columns to numeric
  go_results <- go_results %>%
    mutate(
      Sum_log2FoldChange = as.numeric(Sum_log2FoldChange),
      p.adjust = as.numeric(p.adjust),
      GO_Category = as.factor(GO_Category),  # Convert GO_Category to factor for color mapping
      Wrapped_Description = str_wrap(GO_Description, width = 20)  # Wrap text into multiple lines
    )

  # Handle potential NA values after conversion
  if (any(is.na(go_results$Sum_log2FoldChange) | is.na(go_results$p.adjust))) {
    stop(paste("Error: NA values found in", csv_file))
  }

  # Select the top 10 highest Sum_log2FoldChange values for labeling
  top10 <- go_results %>%
    arrange(desc(Sum_log2FoldChange)) %>%
    slice_head(n = 10)

  # Create the bubble plot (Optimized Label Positioning)
  plot <- ggplot(go_results, aes(x = Sum_log2FoldChange, y = -log10(p.adjust))) +
    geom_point(aes(size = Count, fill = GO_Category), 
               shape = 21,  
               color = "transparent",  
               alpha = 0.8) +  
    geom_text_repel(data = top10, aes(label = Wrapped_Description, color = GO_Category), 
                    size = 3.2, 
                    box.padding = 0.8,  
                    point.padding = 1.2,  
                    force = 3,  
                    nudge_y = 0.2,  
                    nudge_x = 1,  
                    direction = "both",  
                    segment.curvature = -0.3,  
                    max.overlaps = Inf,
                    seed = 42) +  
    scale_fill_manual(values = c("Biological Process" = "#30343F", 
                                 "Molecular Function" = "#656839", 
                                 "Cellular Component" = "#89023E")) +  
    scale_color_manual(values = c("Biological Process" = "#30343F", 
                                  "Molecular Function" = "#656839", 
                                  "Cellular Component" = "#89023E")) +  
    scale_size(range = c(3, 10)) +  
    theme_minimal() +
    theme(
      panel.grid = element_blank(),  
      axis.text.y = element_text(size = 12), 
      axis.text.x = element_text(size = 12),
      legend.key.size = unit(1.5, "cm")  
    ) +
    labs(
      title = paste("GO Enrichment Analysis:", basename(csv_file)),  
      x = "Sum of log2 Fold Change",
      y = "-log10 Adjusted p-value",
      fill = "GO Category",
      color = "GO Category",  
      size = "Gene Count"
    ) +
    guides(
      fill = guide_legend(override.aes = list(size = 6)),  
      color = "none",  
      size = guide_legend(override.aes = list(color = "black"))  
    )

  # Save the plot as an A4-sized PDF
  ggsave(output_pdf, plot, width = 8.27, height = 11.69, units = "in", dpi = 300, device = cairo_pdf)
  
  # Confirmation message
  message("✅ Saved:", output_pdf)
}

# Function to process all CSV files in a directory
process_directory <- function(input_dir, output_dir) {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # Get all CSV files in the directory
  csv_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
  
  # Check if there are any CSV files
  if (length(csv_files) == 0) {
    stop("No CSV files found in the directory:", input_dir)
  }
  
  # Process each CSV file
  for (csv_file in csv_files) {
    # Generate output PDF filename
    pdf_name <- paste0(tools::file_path_sans_ext(basename(csv_file)), "_BubbleChart.pdf")
    output_pdf <- file.path(output_dir, pdf_name)
    
    # Generate the bubble chart
    tryCatch({
      generate_bubble_chart(csv_file, output_pdf)
    }, error = function(e) {
      message("❌ Error processing file:", csv_file, "\n", e$message)
    })
  }
  
  message("🎉 Processing complete! PDFs saved in:", output_dir)
}

# Example Usage:
# Replace `"your_directory_path"` with the actual folder containing CSV files
input_dir <- "GO"
output_dir <- "GO/charts"
process_directory(input_dir, output_dir)
