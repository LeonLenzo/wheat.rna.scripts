# Load Required Libraries
library(clusterProfiler)
library(dplyr)
library(tidyr)
library(stringr)
library(ontologyIndex)

print("✅ Libraries loaded successfully.")

# Load the GO OBO file
go_obo_path <- "GO/go.obo"
if (!file.exists(go_obo_path)) {
  stop("❌ GO OBO file not found. Please download it from: http://purl.obolibrary.org/obo/go.obo")
}

print("📂 Loading GO OBO file...")
go_obo <- get_ontology(go_obo_path, extract_tags = "everything")
print("✅ GO OBO file successfully loaded!")

# Ensure namespace is not NULL
print("🔍 Checking for NULL values in GO namespaces...")
go_obo$namespace[is.null(go_obo$namespace)] <- "Unknown"
print("✅ Namespace check complete.")

# Create GO term mapping
print("📂 Creating GO term mapping...")
go_terms <- data.frame(
  GO_Terms = names(go_obo$id),
  GO_Description = go_obo$name,
  GO_Category = ifelse(is.null(go_obo$namespace), "Unknown", go_obo$namespace),
  stringsAsFactors = FALSE
)
print(paste("✅ GO term mapping created with", nrow(go_terms), "terms."))

# Standardize GO Categories
print("🔍 Standardizing GO Categories...")
go_terms$GO_Category <- case_when(
  go_terms$GO_Category == "biological_process" ~ "Biological Process",
  go_terms$GO_Category == "molecular_function" ~ "Molecular Function",
  go_terms$GO_Category == "cellular_component" ~ "Cellular Component",
  TRUE ~ "Unknown"
)
print("✅ GO Category standardization complete.")

# Load GAF file
gaf_path <- "GO/GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_gene_ontology.gaf"
if (!file.exists(gaf_path)) {
  stop("❌ GO GAF file not found. Please ensure you have downloaded the correct annotation file.")
}

print("📂 Loading GO GAF file... (this may take a moment)")
gaf_raw <- read.table(gaf_path, sep = "\t", quote = "", comment.char = "!", header = FALSE, fill = TRUE, stringsAsFactors = FALSE)
print(paste("✅ GO GAF file successfully loaded with", nrow(gaf_raw), "rows."))

# Define column names based on GAF 2.2 format
print("📝 Assigning column names to GAF file...")
colnames(gaf_raw) <- c("DB", "Gene_ID", "DB_Object_Symbol", "Qualifier", "GO_Terms", "DB_Reference", "Evidence_Code", "With_From", "Aspect", "DB_Object_Name", "DB_Object_Synonym", "DB_Object_Type", "Taxon", "Date", "Assigned_By", "Annotation_Extension", "Gene_Product_Form_ID")

# Extract relevant GO terms and gene mappings
print("📂 Extracting relevant GO terms and gene mappings...")
gaf_clean <- gaf_raw %>%
  select(Gene_ID, GO_Terms, GO_Category = Aspect) %>%
  filter(!is.na(GO_Terms) & GO_Terms != "") %>%
  mutate(
    GO_Category = case_when(
      GO_Category == "F" ~ "Molecular Function",
      GO_Category == "P" ~ "Biological Process",
      GO_Category == "C" ~ "Cellular Component",
      TRUE ~ "Unknown"
    )
  )
print(paste("✅ Extracted", nrow(gaf_clean), "gene-to-GO term mappings."))

# Merge with GO descriptions
print("🔗 Merging GAF data with GO descriptions...")
gaf_clean <- gaf_clean %>%
  left_join(go_terms, by = "GO_Terms")
print(paste("✅ Merged dataset contains", nrow(gaf_clean), "entries."))

# Trim whitespace in all columns
print("🧹 Cleaning whitespace in all columns...")
gaf_clean <- gaf_clean %>%
  mutate(across(everything(), trimws))

# Remove any unwanted or incorrectly duplicated columns
gaf_clean <- gaf_clean %>%
  select(Gene_ID, GO_Terms, GO_Category, GO_Description)  # Keep only relevant columns

# Save cleaned GO annotations
print("💾 Saving cleaned GO annotations...")
write.csv(gaf_clean, "GO_annotations_cleaned.csv", row.names = FALSE)
print("✅ GO annotations cleaned and saved as 'GO_annotations_cleaned.csv'")


# ---- GO Enrichment Analysis ----

# Load cleaned GO annotations
go_clean <- read.csv("GO_XM_annotations_cleaned.csv") # Remember to switch to XM IDs
go_clean$Gene_ID <- trimws(as.character(go_clean$Gene_ID))

go2gene <- go_clean %>% select(GO_Terms, Gene_ID) %>% distinct()
go2desc <- go_clean %>% select(GO_Terms, GO_Description, GO_Category) %>% distinct()

gene_universe <- unique(go_clean$Gene_ID)

# Batch Processing of DESeq2 Files
input_dir <- "deseq/raw_results"
output_dir <- "GO"
if (!dir.exists(output_dir)) dir.create(output_dir)

lfc_threshold <- 2
padj_threshold <- 0.05

deseq_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

for (file in deseq_files) {
  filename <- tools::file_path_sans_ext(basename(file))
  print(paste("🔄 Processing:", filename))

  degs <- read.csv(file)
  degs$gene <- trimws(as.character(degs$gene))
  filtered_genes <- unique(degs$gene[degs$log2FoldChange > lfc_threshold & degs$padj < padj_threshold])
  filtered_genes <- intersect(filtered_genes, gene_universe)
  print(paste("🧬 Filtered genes in", filename, ":", length(filtered_genes)))

  if (length(filtered_genes) == 0) {
    print(paste("⚠️ No DEGs with GO terms found for", filename, "- skipping."))
    next
  }

  ego <- enricher(
    gene = filtered_genes,
    TERM2GENE = go2gene,
    universe = gene_universe,
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1,
    minGSSize = 5
  )

  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    ego_df <- as.data.frame(ego) %>%
      left_join(go2desc, by = c("ID" = "GO_Terms")) %>%
      rename(GO_ID = ID)

    # **Join DEGs with GO Terms to get correct mappings**
    degs_with_go <- degs %>%
      inner_join(go2gene, by = c("gene" = "Gene_ID"))

    # **Summarize log2FoldChange per GO term**
    lfc_sums <- degs_with_go %>%
      group_by(GO_Terms) %>%
      summarise(Sum_log2FoldChange = sum(log2FoldChange, na.rm = TRUE), .groups = "drop")

    ego_df <- ego_df %>%
      left_join(lfc_sums, by = c("GO_ID" = "GO_Terms"))

    output_file <- file.path(output_dir, paste0(filename, "_GO_enrichment.csv"))
    write.csv(ego_df, output_file, row.names = FALSE)
    
    print(paste("✅ Enrichment results saved to:", output_file))
  } else {
    print(paste("⚠️ No significant GO terms found for", filename))
  }
}

print("🎯 Batch processing complete! All results saved.")
