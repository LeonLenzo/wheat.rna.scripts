# Load necessary libraries
library(dplyr)
library(ggplot2)
library(tidyr)

# Ensure required directories exist
if (!dir.exists("GO")) {
  dir.create("GO")
}

# Load DESeq2 results
deseq_file <- "deseq/raw_results/deseq2_SN15_vs_Tween.csv"
if (!file.exists(deseq_file)) {
  stop(paste("ERROR: File not found:", deseq_file))
}
deseq_results <- read.csv(deseq_file)
colnames(deseq_results) <- tolower(colnames(deseq_results))

# Load custom GO annotations
go_file <- "mapping/GO_cds_annotations_cleaned.csv"
if (!file.exists(go_file)) {
  stop(paste("ERROR: File not found:", go_file))
}
custom_go <- read.csv(go_file)

# Filter genes based on fold change thresholds
upregulated <- deseq_results %>%
  filter(log2foldchange > 3 & padj < 0.01) %>%
  pull(gene)  

downregulated <- deseq_results %>%
  filter(log2foldchange < -3 & padj < 0.01) %>%
  pull(gene)

# Merge custom GO annotations with DESeq2 gene lists
up_go <- custom_go %>% filter(Gene_ID %in% upregulated)
down_go <- custom_go %>% filter(Gene_ID %in% downregulated)

# Count occurrences of GO terms
go_up_count <- up_go %>%
  group_by(GO_Terms, GO_Description, GO_Category) %>%
  summarise(up_count = n(), .groups = 'drop')

go_down_count <- down_go %>%
  group_by(GO_Terms, GO_Description, GO_Category) %>%
  summarise(down_count = n(), .groups = 'drop')

# Get total GO term background
all_go_count <- custom_go %>%
  group_by(GO_Terms) %>%
  summarise(total_count = n(), .groups = 'drop')

# Fisher's exact test function
fisher_test <- function(go_count, all_go_count, gene_count) {
  if (nrow(go_count) == 0) return(data.frame(GO_Terms = character(), GO_Description = character(), GO_Category = character(), adj_p_value = numeric()))

  enriched <- go_count %>%
    left_join(all_go_count, by = "GO_Terms") %>%
    rowwise() %>%
    mutate(
      p_value = tryCatch(
        fisher.test(matrix(c(n(), max(1, gene_count - n()),
                             total_count - n(), max(1, sum(all_go_count$total_count) - total_count)),
                           nrow = 2))$p.value,
        error = function(e) NA
      )
    ) %>%
    ungroup() %>%
    mutate(adj_p_value = ifelse(is.na(p_value), NA, p.adjust(p_value, method = "BH"))) %>%
    filter(!is.na(p_value) & adj_p_value < 0.05) %>%
    arrange(adj_p_value)

  return(enriched)
}

# Apply Fisher's test separately, ensuring non-empty dataframes
go_up_enrichment <- fisher_test(go_up_count, all_go_count, length(upregulated))
go_down_enrichment <- fisher_test(go_down_count, all_go_count, length(downregulated))

if (nrow(go_up_enrichment) == 0) {
  go_up_enrichment <- data.frame(GO_Terms = character(), GO_Description = character(), GO_Category = character(), adj_p_value = numeric())
}
if (nrow(go_down_enrichment) == 0) {
  go_down_enrichment <- data.frame(GO_Terms = character(), GO_Description = character(), GO_Category = character(), adj_p_value = numeric())
}

# Merge upregulated and downregulated enrichment results into a single table
combined_go_enrichment <- full_join(
  go_up_enrichment %>% rename(up_adj_p_value = adj_p_value),
  go_down_enrichment %>% rename(down_adj_p_value = adj_p_value),
  by = c("GO_Terms", "GO_Description", "GO_Category")
) %>%
  mutate(up_adj_p_value = replace_na(up_adj_p_value, 1),
         down_adj_p_value = replace_na(down_adj_p_value, 1)) %>%
  arrange(up_adj_p_value, down_adj_p_value)

# Save results
if (nrow(go_up_enrichment) > 0) {
  write.csv(go_up_enrichment, "GO/GO_enrichment_upregulated.csv", row.names = FALSE)
}
if (nrow(go_down_enrichment) > 0) {
  write.csv(go_down_enrichment, "GO/GO_enrichment_downregulated.csv", row.names = FALSE)
}
if (nrow(combined_go_enrichment) > 0) {
  write.csv(combined_go_enrichment, "GO/GO_enrichment_combined.csv", row.names = FALSE)
}

# Create summary table for bar plot: Total up/down per GO Category
category_counts <- full_join(
  go_up_count %>% group_by(GO_Category) %>% summarise(up_count = sum(up_count), .groups = "drop"),
  go_down_count %>% group_by(GO_Category) %>% summarise(down_count = sum(down_count), .groups = "drop"),
  by = "GO_Category"
) %>%
  mutate(up_count = replace_na(up_count, 0),
         down_count = replace_na(down_count, 0))

# Save category summary table
if (nrow(category_counts) > 0) {
  write.csv(category_counts, "GO/GO_category_summary.csv", row.names = FALSE)
}

# Visualization: Save plots as PNG instead of PDF
if (nrow(combined_go_enrichment) > 0) {
  png("GO/GO_Enrichment.png", width = 1200, height = 800, res = 150)
  print(
    ggplot(combined_go_enrichment[1:min(20, nrow(combined_go_enrichment)), ],
           aes(x = reorder(GO_Description, -pmin(up_adj_p_value, down_adj_p_value)), 
               y = -log10(pmin(up_adj_p_value, down_adj_p_value)), fill = GO_Category)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      labs(title = "GO Enrichment (Combined Up & Down)", x = "GO Term", y = "-log10(Adjusted p-value)") +
      theme_minimal()
  )
  dev.off()
}

# New Bar Plot: Total Up/Down Per GO Category
if (nrow(category_counts) > 0) {
  png("GO/GO_Category_Counts.png", width = 1200, height = 800, res = 150)
  print(
    ggplot(category_counts, aes(x = GO_Category)) +
      geom_bar(aes(y = up_count, fill = "Upregulated"), stat = "identity", position = "dodge") +
      geom_bar(aes(y = -down_count, fill = "Downregulated"), stat = "identity", position = "dodge") +
      scale_y_continuous(labels = abs) +  # Show absolute values on y-axis
      labs(title = "Total Up & Down-Regulated Genes Per GO Category", x = "GO Category", y = "Gene Count") +
      scale_fill_manual(values = c("Upregulated" = "blue", "Downregulated" = "red")) +
      theme_minimal() +
      coord_flip()
  )
  dev.off()
}
