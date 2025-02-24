### Step 1: Load Libraries ###
library(clusterProfiler)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)

### Step 2: Load GO Annotations (Static File) ###
go_clean <- read_csv("mapping/GO2_annotations_cleaned.csv")  # Update with actual filename

# Ensure Gene ID formatting consistency
go_clean <- go_clean %>%
  mutate(
    Gene_ID = trimws(as.character(Gene_ID)),
    GO_Terms = as.character(GO_Terms),
    GO_Description = trimws(GO_Description)  # Remove leading/trailing spaces
  )

# Define the GO annotation mapping (TERM2GENE)
go2gene <- go_clean %>%
  select(GO_Terms, Gene_ID) %>%
  distinct()

# Create a GO Term to Description mapping
go2desc <- go_clean %>%
  select(GO_Terms, GO_Description, GO_Category) %>%
  distinct()

# Define the gene universe (all genes with GO terms)
gene_universe <- unique(go_clean$Gene_ID)

### Step 3: Batch Processing of DESeq2 Files ###
input_dir <- "deseq/filtered_results"  # Set directory containing DESeq2 outputs
output_dir <- "GO/GO2"  # Set directory for enrichment results

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Define filtering thresholds (adjustable)
lfc_threshold <- 10  # Log2 Fold Change threshold
padj_threshold <- 0.01  # Adjusted p-value threshold

# List all CSV files in the input directory
deseq_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

# Loop through each DESeq2 file
for (file in deseq_files) {
  
  # Extract filename for labeling output
  filename <- tools::file_path_sans_ext(basename(file))
  
  print(paste("🔄 Processing:", filename))
  
  ### Step 4: Load DESeq2 Results ###
  degs <- read_csv(file)
  
  print(paste("📄 Loaded:", filename, "with", nrow(degs), "rows"))
  
  # Ensure Gene ID formats are consistent
  degs <- degs %>%
    mutate(gene = trimws(as.character(gene)))
  
  # **Filter genes dynamically based on threshold values**
  filtered_genes <- degs %>%
    filter(log2FoldChange > lfc_threshold & padj < padj_threshold) %>%
    pull(gene) %>%
    unique()
  
  # Ensure we're only keeping genes present in the GO annotation file
  filtered_genes <- intersect(filtered_genes, gene_universe)
  
  print(paste("🧬 Filtered genes in", filename, ":", length(filtered_genes)))

  if (length(filtered_genes) == 0) {
    print(paste("⚠️ No DEGs with GO terms found for", filename, "- skipping."))
    next
  }

  ### Step 5: Perform GO Enrichment ###
  ego <- enricher(
    gene = filtered_genes,  
    TERM2GENE = go2gene, 
    universe = gene_universe,  
    pvalueCutoff = 0.05,  
    qvalueCutoff = 0.1,  
    minGSSize = 5  
  )

  # Check if enrichment was successful
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    
    # Convert enrichment results to a dataframe
    ego_df <- as.data.frame(ego)

    # Ensure ID column is character before merging
    ego_df <- ego_df %>%
      mutate(ID = as.character(ID)) %>%
      left_join(go2desc, by = c("ID" = "GO_Terms")) %>%
      rename(GO_ID = ID)

    ### Step 6: Compute Sum of log2FoldChange for Each GO Term ###
    # Merge DESeq2 results with GO annotations
    degs_with_go <- degs %>%
      inner_join(go_clean %>% select(Gene_ID, GO_Terms), by = c("gene" = "Gene_ID"))

    # Sum log2FoldChange for each GO term
    lfc_sums <- degs_with_go %>%
      group_by(GO_Terms) %>%
      summarise(Sum_log2FoldChange = sum(log2FoldChange, na.rm = TRUE), .groups = "drop")

    # Merge log2FoldChange sums into enrichment results
    ego_df <- ego_df %>%
      left_join(lfc_sums, by = c("GO_ID" = "GO_Terms"))

    # Save results
    output_file <- file.path(output_dir, paste0(filename, "_GO_enrichment.csv"))
    write_csv(ego_df, output_file)
    
    print(paste("✅ Enrichment results saved to:", output_file))
  } else {
    print(paste("⚠️ No significant GO terms found for", filename))
  }
}

print("🎯 Batch processing complete! All results saved.")
