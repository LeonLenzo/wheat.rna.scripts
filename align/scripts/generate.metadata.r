# Load required libraries
library(dplyr)

# Define the directory containing Kallisto output
kallisto_dir <- "./quant/cs_cds" 

# List all abundance.h5 files recursively in the directory
files <- list.files(kallisto_dir, pattern = "abundance.tsv", recursive = TRUE, full.names = TRUE)

# Extract sample names from the directory structure or file names
samples <- basename(dirname(files))

# Define a mapping of sample names to conditions
# Modify this list as per your experimental design
condition_mapping <- data.frame(
  sample = samples,
  condition = case_when(
    grepl("Tween", samples, ignore.case = TRUE) ~ "Tween",
    grepl("SN15", samples, ignore.case = TRUE) ~ "SN15",
    grepl("Pf2", samples, ignore.case = TRUE) ~ "Pf2",
    #grepl("3KO", samples, ignore.case = TRUE) ~ "3KO",
    TRUE ~ "unknown" # Assign "unknown" for unmatched samples
  )
)

# Combine sample names, conditions, and paths into a metadata data frame
metadata <- data.frame(
  sample = samples,
  condition = condition_mapping$condition,
  path = files,
  stringsAsFactors = FALSE
)

# View the metadata table
print(metadata)

# Save the metadata to a CSV file
output_path <- file.path(kallisto_dir, "metadata-no3ko.csv")
write.csv(metadata, output_path, row.names = FALSE)

cat("Metadata file created at:", output_path, "\n")
