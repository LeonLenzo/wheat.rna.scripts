# Load Required Libraries
library(dplyr)
library(tidyr)
library(stringr)

# Load the GO annotation file (Tab-Separated File)
go_raw <- read.delim("/home/leonl/wheat.rna/mapping/IWGSC/GO.genes.txt", header = FALSE, stringsAsFactors = FALSE)

# Rename columns
colnames(go_raw) <- c("prot_ID", "Description", "GO_Terms")

# Function to extract and organize GO terms into separate categories
extract_go_terms <- function(data) {
  data %>%
    filter(!is.na(GO_Terms) & GO_Terms != "") %>%  # Remove empty GO term values
    separate_rows(GO_Terms, sep = ";") %>%  # Split multiple GO terms per gene
    mutate(
      GO_ID = str_extract(GO_Terms, "GO:\\d+"),  # Extract GO ID
      GO_Description = str_remove(GO_Terms, "GO:\\d+ "),  # Extract GO description
      GO_Category = case_when(
        str_detect(GO_Terms, "MF:") ~ "Molecular Function",
        str_detect(GO_Terms, "BP:") ~ "Biological Process",
        str_detect(GO_Terms, "CC:") ~ "Cellular Component",
        TRUE ~ "Unknown"
      )
    ) %>%
    filter(!is.na(GO_ID))  # Remove rows without valid GO terms
}

# Apply function to clean GO data
go_clean <- extract_go_terms(go_raw)

# Reshape data to wide format: One row per gene with GO categories as separate columns
go_shiny <- go_clean %>%
  group_by(prot_ID, Description, GO_Category) %>%
  summarise(GO_Terms = paste0(GO_Description, " [", GO_ID, "]", collapse = "; "), .groups = "drop") %>%
  pivot_wider(names_from = GO_Category, values_from = GO_Terms, values_fill = "")

# Ensure columns are in the correct order
go_shiny <- go_shiny %>%
  select(prot_ID, Description, `Molecular Function`, `Biological Process`, `Cellular Component`)

# Save as a CSV file
write.csv(go_shiny, "mapping/GO_annotations_shiny.csv", row.names = FALSE)

print("✅ GO annotations reformatted for Shiny app and saved as 'GO_annotations_shiny.csv'.")
