# Load necessary libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(plotly)

# Ensure required directories exist
if (!dir.exists("GO")) {
  dir.create("GO")
}

# Define the directory containing DESeq2 result files
deseq_dir <- "deseq/raw_results/"

# Get all CSV files in the directory
deseq_files <- list.files(path = deseq_dir, pattern = "*.csv", full.names = TRUE)

# Load custom GO annotations
go_file <- "mapping/GO_cds_annotations_cleaned.csv"
if (!file.exists(go_file)) {
  stop(paste("ERROR: File not found:", go_file))
}
custom_go <- read.csv(go_file)

# Initialize a list to store GO enrichment results
all_conditions_data <- list()

# Process each DESeq dataset
for (file in deseq_files) {
  condition_label <- tools::file_path_sans_ext(basename(file))  # Extract filename without extension

  deseq_results <- read.csv(file)
  colnames(deseq_results) <- tolower(colnames(deseq_results))

  # Filter genes based on fold change thresholds
  upregulated <- deseq_results %>%
    filter(log2foldchange > 2 & padj < 0.01) %>%
    pull(gene)

  downregulated <- deseq_results %>%
    filter(log2foldchange < -2 & padj < 0.01) %>%
    pull(gene)

  # Merge GO annotations
  up_go <- custom_go %>% filter(Gene_ID %in% upregulated)
  down_go <- custom_go %>% filter(Gene_ID %in% downregulated)

  # Count occurrences of GO terms
  go_up_count <- up_go %>%
    group_by(GO_Terms, GO_Description, GO_Category) %>%
    summarise(up_count = n(), .groups = 'drop')

  go_down_count <- down_go %>%
    group_by(GO_Terms, GO_Description, GO_Category) %>%
    summarise(down_count = n(), .groups = 'drop')

  # Merge up/down enrichment results (keeping counts)
  combined_go_enrichment <- full_join(
    go_up_count %>% rename(Upregulated = up_count),
    go_down_count %>% rename(Downregulated = down_count),
    by = c("GO_Terms", "GO_Description", "GO_Category")
  ) %>%
    mutate(
      Upregulated = replace_na(Upregulated, 0),
      Downregulated = replace_na(Downregulated, 0),
      Condition = condition_label  # Add condition column
    )

  all_conditions_data[[condition_label]] <- combined_go_enrichment
}

# Combine all conditions into a single table
final_combined_go_enrichment <- bind_rows(all_conditions_data)

# Save final combined enrichment table
write.csv(final_combined_go_enrichment, "GO/GO_enrichment_combined_multiple_conditions.csv", row.names = FALSE)

# Generate separate plots for each GO category, comparing conditions
unique_categories <- unique(final_combined_go_enrichment$GO_Category)

for (cat in unique_categories) {
  category_data <- final_combined_go_enrichment %>%
    filter(GO_Category == cat) %>%
    pivot_longer(cols = c("Upregulated", "Downregulated"), 
                 names_to = "Regulation", values_to = "Gene_Count") %>%
    mutate(Regulation = factor(Regulation, levels = c("Downregulated", "Upregulated")))

  if (nrow(category_data) > 0) {
    # Adjust height with a max cap to prevent Cairo crash
    height_factor <- min(max(1200, nrow(category_data) * 30), 8000)

    png(paste0("GO/GO_", gsub(" ", "_", cat), "_Multiple_Conditions.png"), 
        width = 2500, height = height_factor, res = 200)

    print(
      ggplot(category_data, aes(x = reorder(GO_Description, Gene_Count), y = Gene_Count, fill = Regulation)) +
        geom_bar(stat = "identity", position = "dodge") +  # Dodge to separate up/down bars
        facet_grid(. ~ Condition, scales = "free_x", space = "free_x") +  # Side-by-side condition comparison
        coord_flip() +
        scale_fill_manual(values = c("Upregulated" = "blue", "Downregulated" = "red")) +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 40)) +  # Wrap long GO term labels
        labs(title = paste("GO Enrichment -", cat), 
             x = "GO Term", 
             y = "Number of Genes") +
        theme_minimal(base_size = ifelse(nrow(category_data) > 200, 10, 14))  # Reduce font size for large plots
    )

    dev.off()

    
  }
}

for (cat in unique_categories) {
  category_data <- final_combined_go_enrichment %>%
    filter(GO_Category == cat) %>%
    mutate(Total_Genes = Upregulated + Downregulated) %>%  # Calculate total gene count
    arrange(desc(Total_Genes)) %>%  # Sort GO terms by total gene count
    slice_head(n = 50) %>%  # Keep only the top 50 GO terms per category
    pivot_longer(cols = c("Upregulated", "Downregulated"), 
                 names_to = "Regulation", values_to = "Gene_Count") %>%
    mutate(Regulation = factor(Regulation, levels = c("Downregulated", "Upregulated")))

  if (nrow(category_data) > 0) {
    # Create interactive ggplot
    p <- ggplot(category_data, aes(x = reorder(GO_Description, Gene_Count), y = Gene_Count, fill = Regulation, text = paste("GO Term:", GO_Description, "<br>Genes:", Gene_Count))) +
      geom_bar(stat = "identity", position = "dodge") +
      facet_grid(. ~ Condition, scales = "free_x", space = "free_x") +
      coord_flip() +
      scale_fill_manual(values = c("Upregulated" = "#990000", "Downregulated" = "#156082")) +
      scale_x_discrete(labels = function(x) str_wrap(x, width = 40)) +
      labs(title = paste("GO Enrichment -", cat, "(Top 50)"), 
           x = "GO Term", 
           y = "Number of Genes") +
      theme_minimal(base_size = ifelse(nrow(category_data) > 200, 10, 14))

    # Convert ggplot to an interactive plotly chart
    interactive_plot <- ggplotly(p, tooltip = "text")  # Enable hover text

    # Save the interactive plot as an HTML file
    htmlwidgets::saveWidget(interactive_plot, file = paste0("GO/GO_", gsub(" ", "_", cat), "_Top50_Multiple_Conditions.html"))
  }
}