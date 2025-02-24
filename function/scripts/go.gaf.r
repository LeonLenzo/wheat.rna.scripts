library(dplyr)
library(tidyr)
library(stringr)
library(readr)

# Load the GAF file (modify filename as needed)
gaf_file <- "GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_gene_ontology.gaf"

# Read the .gaf file, skipping the header lines (which start with '!')
go_annotations <- read_tsv(
  gaf_file, comment = "!", col_names = FALSE, show_col_types = FALSE
)

# View the first few rows
print(head(go_annotations))

go_clean <- go_annotations %>%
  select(Gene_ID = X2, GO_Terms = X5, GO_Category = X9) %>%
  mutate(
    GO_Category = case_when(
      GO_Category == "F" ~ "Molecular Function",
      GO_Category == "P" ~ "Biological Process",
      GO_Category == "C" ~ "Cellular Component",
      TRUE ~ GO_Category
    )
  ) %>%
  distinct()

# Display cleaned GO mappings
print(head(go_clean, 10))

# Save cleaned annotations
write.csv(go_clean, "GO_annotations_cleaned.csv", row.names = FALSE)

print("✅ GO annotations cleaned and saved as 'GO_annotations_cleaned.csv'")

# Load the cleaned GO annotation file
go_clean <- read.csv("GO_annotations_cleaned.csv")

# Ensure Gene ID formatting consistency
go_clean$Gene_ID <- trimws(as.character(go_clean$Gene_ID))

# Define GO annotation mappings
go2gene <- go_clean %>% 
  select(GO_Terms, Gene_ID) %>% 
  distinct()

# Create GO Term to Description mapping
go2desc <- go_clean %>%
  select(GO_Terms, GO_Category) %>%
  distinct()

