# Load Required Libraries
library(dplyr)
library(tidyr)
library(stringr)

# Load the GO annotation file
go_raw <- read.csv("mapping/GO_cds.csv", stringsAsFactors = FALSE)  # Ensure no factors

# Check column names
print("📄 Available Columns in Input File:")
print(colnames(go_raw))

# Ensure expected columns exist
expected_columns <- c("Transcript", "MF", "BP", "CC")
missing_cols <- setdiff(expected_columns, colnames(go_raw))

if (length(missing_cols) > 0) {
  stop(paste("❌ Missing required columns:", paste(missing_cols, collapse = ", ")))
}

# Function to clean and extract GO terms from a given column (MF, BP, or CC)
extract_go_terms <- function(data, column_name, category_label) {
  data %>%
    select(Gene_ID = Transcript, GO_Terms = !!sym(column_name)) %>%
    filter(!is.na(GO_Terms) & GO_Terms != "") %>%  # Remove empty values
    separate_rows(GO_Terms, sep = ";") %>%  # Split multiple GO terms per gene
    mutate(
      GO_Description = str_extract(GO_Terms, "^[^\\[]+"),  # Extract text before "[GO:"
      GO_Terms = str_extract(GO_Terms, "(?<=\\[GO:)\\d+(?=\\])"),  # Extract GO ID inside "[GO: ]"
      GO_Category = category_label  # Assign category (MF, BP, or CC)
    ) %>%
    select(Gene_ID, GO_Terms, GO_Description, GO_Category) %>%
    filter(!is.na(GO_Terms) & GO_Terms != "")  # Ensure valid GO terms
}

# Process each GO category separately
go_mf <- extract_go_terms(go_raw, "MF", "Molecular Function")
go_bp <- extract_go_terms(go_raw, "BP", "Biological Process")
go_cc <- extract_go_terms(go_raw, "CC", "Cellular Component")

# Combine all cleaned GO terms
go_clean <- bind_rows(go_mf, go_bp, go_cc) %>%
  distinct()  # Ensure uniqueness

# Trim whitespace in all columns
go_clean <- go_clean %>%
  mutate(across(everything(), trimws))

# Debugging: Check structure
print("🔍 Sample of cleaned GO annotation data:")
print(head(go_clean, 10))

# Check for duplicate GO entries
print("🔍 Checking for duplicated GO terms:")
print(sum(duplicated(go_clean$GO_Terms)))

# Save cleaned file
write.csv(go_clean, "mapping/GO_cds_annotations_cleaned.csv", row.names = FALSE)

print("✅ GO annotations cleaned and saved as 'GO_annotations_cleaned.csv'")
