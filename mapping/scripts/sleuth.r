# Load necessary libraries
library(sleuth)

# Step 1: Load Metadata
# Ensure your metadata file has columns: sample, condition, path
metadata <- read.csv("metadata.csv", stringsAsFactors = FALSE)

# Inspect metadata
print(metadata)

# Step 2: Prepare Sleuth Object
# Lower the filtering threshold if needed (default min_reads = 5)
so <- sleuth_prep(metadata, ~condition, min_reads = 5)

# Step 3: Fit Models
# Full model includes the condition; reduced model assumes no effect of condition
so <- sleuth_fit(so, ~condition, "full")
so <- sleuth_fit(so, ~1, "reduced")

# Step 4: Likelihood Ratio Test (LRT)
# Compare full and reduced models
so <- sleuth_lrt(so, "reduced", "full")

# Step 5: Retrieve Results
# Extract differential expression results
results_table <- sleuth_results(so, "reduced:full", "lrt")

# Save results to a CSV file
write.csv(results_table, "sleuth_lrt_results.csv", row.names = FALSE)

# Step 6: Filter Significant Results
# Filter transcripts with q-value < 0.05 (significant after FDR correction)
significant <- results_table[results_table$qval < 0.05, ]

# Save significant results to a CSV file
write.csv(significant, "significant_differential_expression.csv", row.names = FALSE)

# Step 7: PCA Plot
# Visualize sample clustering using PCA
pca_plot <- plot_pca(so, color_by = "condition", return_data = TRUE)
ggsave("pca_plot.png", plot = pca_plot, width = 8, height = 6, dpi = 300)

# Step 8: Examine Individual Transcripts
# Visualize expression of a specific transcript (replace "target_id" with your gene of interest)
plot_bootstrap(so, target_id = "TraesMAC1A03G00000010.1", units = "tpm", color_by = "condition")

# Step 9: Wald Test (Optional)
# If interested in pairwise comparisons, perform a Wald test
so <- sleuth_wt(so, "conditionSN15")  # Replace "Tween" with your condition of interest
wald_results <- sleuth_results(so, "conditionSN15", "wt")

# Save Wald test results to a CSV file
write.csv(wald_results, "wald_test_results.csv", row.names = FALSE)

# Step 10: Interactive Exploration
# Use Sleuth's web-based interface to explore results interactively
sleuth_live(so)

# Step 11: Add Directionality (Optional)
# Calculate directionality of differential expression (log2 fold changes)
# Extract TPM values from Sleuth object
tpm_table <- sleuth_to_matrix(so, "obs_norm", "tpm")

# Calculate mean TPM for each condition
mean_condition1 <- rowMeans(tpm_table[, metadata$condition == "Tween"])
mean_condition2 <- rowMeans(tpm_table[, metadata$condition == "SN15"])
mean_condition3 <- rowMeans(tpm_table[, metadata$condition == "Pf2"])

# Combine means into a data frame
mean_tpm <- data.frame(
  target_id = rownames(tpm_table),
  mean_condition1 = mean_condition1,
  mean_condition2 = mean_condition2,
  mean_condition3 = mean_condition3
)

# Determine directionality
mean_tpm$direction <- apply(mean_tpm[, -1], 1, function(x) {
  conditions <- c("Tween", "SN15", "Pf2")
  conditions[which.max(x)]  # Identify the condition with the highest mean TPM
})

# Add directionality to results table
results_table <- merge(results_table, mean_tpm[, c("target_id", "direction")], by = "target_id")

# Save table with directionality to CSV
write.csv(results_table, "sleuth_with_directionality.csv", row.names = FALSE)


# Save table with directionality to CSV
write.csv(results_table, "sleuth_with_log_fc.csv", row.names = FALSE)
