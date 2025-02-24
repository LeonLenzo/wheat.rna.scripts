# Ensure biomaRt is installed and loaded
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("biomaRt")
library(biomaRt)

# ---------------- Step 1: Connect to Ensembl Plants ----------------
ensembl_plants <- useEnsemblGenomes(biomart = "plants_mart")

# Load wheat dataset (Triticum aestivum)
mart <- useDataset("taestivum_eg_gene", mart = ensembl_plants)

cat("Successfully connected to Ensembl Plants BioMart.\n")

# ---------------- Step 2: Load Attribute List from CSV ----------------
attribute_file <- "mapping/selected_attributes.csv"  # Update if needed

if (!file.exists(attribute_file)) {
  stop(paste("Error: Attribute file", attribute_file, "not found."))
}

# Read the attribute names from CSV (assuming it's a single column without a header)
selected_attributes <- read.csv(attribute_file, header = FALSE, stringsAsFactors = FALSE)[,1]

# Ensure attributes are loaded
if (length(selected_attributes) == 0) {
  stop("Error: No attributes found in the selected_attributes.csv file.")
}

cat("Fetching", length(selected_attributes), "attributes one by one...\n")

# ---------------- Step 3: Retrieve Data One Attribute at a Time ----------------
temp_files <- c()  # Store file paths for later merging

for (attribute in selected_attributes) {
  cat("Fetching:", attribute, "\n")
  
  # Fetch data for the current attribute
  result <- tryCatch({
    getBM(attributes = c("ensembl_gene_id", attribute), mart = mart)
  }, error = function(e) {
    cat("Error fetching", attribute, "- skipping.\n")
    return(NULL)
  })
  
  # Skip if no results
  if (is.null(result) || nrow(result) == 0) next

  # Save each attribute as a temporary file
  temp_file <- paste0("temp_", attribute, ".csv")
  write.csv(result, temp_file, row.names = FALSE)
  temp_files <- c(temp_files, temp_file)
}

# ---------------- Step 4: Merge All Data ----------------
cat("Merging all retrieved attributes...\n")

# Read and merge all saved files
merged_data <- Reduce(function(x, y) merge(x, read.csv(y), by = "ensembl_gene_id", all = TRUE),
                      lapply(temp_files, read.csv))

# ---------------- Step 5: Save Final Merged Data ----------------
output_file <- "wheat_selected_attributes_combined.csv"
write.csv(merged_data, output_file, row.names = FALSE)

cat("Completed! Merged wheat annotations saved to", output_file, "\n")

# ---------------- Step 6: Clean Up Temporary Files ----------------
file.remove(temp_files)
cat("Temporary files deleted.\n")
