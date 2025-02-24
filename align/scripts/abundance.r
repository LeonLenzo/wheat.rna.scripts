# Load required library
library(dplyr)

# Define the directory containing the abundance files
abundance_dir <- "./quant"

# Get a list of abundance.tsv files
abundance_files <- list.files(abundance_dir, pattern = "abundance.tsv", full.names = TRUE)

# Initialize an empty data frame to store combined data
combined_table <- NULL

# Iterate through each file and extract the required column
for (file in abundance_files) {
  # Read the abundance.tsv file
  abundance_data <- read.delim(file, header = TRUE, stringsAsFactors = FALSE)
  
  # Extract the sample ID from the file path
  sample_id <- basename(dirname(file))
  
  # Select target_id (transcript ID) and the desired metric (e.g., est_counts or tpm)
  sample_table <- abundance_data %>%
    select(target_id, est_counts) %>%
    rename(!!sample_id := est_counts)  # Rename the est_counts column to the sample ID
  
  # Join the current sample's data with the combined table
  if (is.null(combined_table)) {
    combined_table <- sample_table
  } else {
    combined_table <- full_join(combined_table, sample_table, by = "target_id")
  }
}

# Save the combined table to a CSV file
write.csv(combined_table, "combined_abundance_table.csv", row.names = FALSE)

# Print the first few rows of the combined table
head(combined_table)
