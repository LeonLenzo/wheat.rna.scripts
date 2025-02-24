library(biomaRt)

# ---------------- Step 1: Load RefSeq mRNA Predicted IDs from CSV ----------------
# Define file paths
input_file <- "mapping/refseq_mrna_predicted_id.csv"  # Ensure this file exists
output_file <- "mapping/refseq_to_ensembl_mapping.csv"

# Check if file exists
if (!file.exists(input_file)) {
  stop(paste("Error: Input file", input_file, "not found. Please check the path."))
}

# Read the CSV file (assuming it has a column "refseq_mrna_predicted_id")
refseq_data <- read.csv(input_file, stringsAsFactors = FALSE)

# Ensure the column exists
if (!"refseq_mrna_predicted_id" %in% colnames(refseq_data)) {
  stop("Error: The CSV file must have a column named 'refseq_mrna_predicted_id'.")
}

# Extract unique RefSeq mRNA predicted IDs
refseq_ids <- unique(na.omit(refseq_data$refseq_mrna_predicted_id))  # Remove duplicates & NAs

# Check if the list is empty
if (length(refseq_ids) == 0) {
  stop("Error: No valid RefSeq mRNA Predicted IDs found in the input file.")
}

cat("Loaded", length(refseq_ids), "RefSeq mRNA Predicted IDs.\n")

# ---------------- Step 2: Connect to Ensembl Plants BioMart ----------------
ensembl_plants <- useEnsemblGenomes(biomart = "plants_mart")

# Load wheat dataset (Triticum aestivum)
mart <- useDataset("taestivum_eg_gene", mart = ensembl_plants)

cat("Successfully connected to Ensembl Plants BioMart.\n")

# ---------------- Step 3: Retrieve Ensembl Gene IDs ----------------
cat("Querying Ensembl for gene mappings...\n")

# Query Ensembl for mappings using the correct filter
gene_mapping <- getBM(
  attributes = c("refseq_mrna_predicted", "ensembl_gene_id"),
  filters = "refseq_mrna_predicted",
  values = refseq_ids,
  mart = mart
)

# ---------------- Step 4: Save Results to CSV ----------------
# Merge with original data for better tracking
merged_data <- merge(refseq_data, gene_mapping, by.x = "refseq_mrna_predicted_id", by.y = "refseq_mrna_predicted", all.x = TRUE)

# Save to CSV
write.csv(merged_data, output_file, row.names = FALSE)

cat("Mapping completed! Results saved to", output_file, "\n")
