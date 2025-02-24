# Load required libraries
library(DESeq2)
library(tximport)
library(readr)
library(dplyr)  
library(fs)  
library(ggplot2)

# Load metadata
metadata <- read.csv("quant/cs_cds/metadata.csv")
metadata$condition <- factor(metadata$condition)

# Create named vector of file paths
files <- metadata$path
names(files) <- metadata$sample

# Import Kallisto quantifications
txi <- tximport(files, type = "kallisto", txOut = TRUE)

# Create DESeq2 dataset
dds <- DESeqDataSetFromTximport(txi, colData = metadata, design = ~ condition)

# Check total genes before filtering
nrow(dds)
# Apply filtering
dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]
# Check total genes after filtering
nrow(dds)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Define filtering thresholds
log2fc_threshold <- 2  
high_log2fc_threshold <- 10  
padj_threshold <- 0.05
high_padj_threshold <- 1e-10  # More relaxed padj for high log2FC genes

# Ensure output directories exist
dir_create("deseq/raw_results")
dir_create("deseq/filtered_results")
dir_create("deseq/high_l2fc_filtered_results")
dir_create("deseq/plots")

# PCA Plot for All Samples
vsd <- vst(dds, blind = TRUE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = 4, shape = 21, fill = "white") +
  geom_text(vjust = -1, size = 3) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  ggtitle("PCA Plot of Samples") +
  theme_minimal()

ggsave("deseq/plots/PCA_plot.png", plot = p, width = 8, height = 6, dpi = 300)

# Get all condition levels
conditions <- levels(metadata$condition)

# 🔹 **Loop through all pairwise comparisons (both orders)**
for (i in seq_along(conditions)) {
  for (j in seq_along(conditions)) {
    if (i != j) {  # Allow both (A vs B) and (B vs A)
      numerator <- conditions[j]
      denominator <- conditions[i]
      comparison_name <- paste(numerator, "_vs_", denominator, sep = "")

      # Extract all DESeq2 results (keep as DESeqResults object)
      res_deseq <- results(dds, contrast = c("condition", numerator, denominator))

      # Convert to data frame for processing
      res <- as.data.frame(res_deseq)
      res$gene <- rownames(res)

      # Extract clean gene names (remove lcl|NC_... prefixes)
      res$gene <- gsub("^.*_cds_(XP_\\d+).*", "\\1", res$gene)

      # Remove non-coding transcripts (XR-prefixed)
      res <- res %>% filter(!grepl("^XR", gene))

      # Save full DESeq2 output
      write.csv(res, file = paste0("deseq/raw_results/deseq2_", comparison_name, ".csv"), row.names = FALSE)

      # Apply standard filtering (log2FC & padj)
      res_filtered <- res %>% filter(abs(log2FoldChange) >= log2fc_threshold & padj <= padj_threshold)
      write.csv(res_filtered, file = paste0("deseq/filtered_results/deseq2_", comparison_name, "_filtered.csv"), row.names = FALSE)

      # 🔹 **Filter for high |log2FC| genes with a higher padj cutoff**
      res_high_l2fc <- res %>% filter(abs(log2FoldChange) >= high_log2fc_threshold & padj <= high_padj_threshold)
      write.csv(res_high_l2fc, file = paste0("deseq/high_l2fc_filtered_results/deseq2_", comparison_name, "_high_l2fc.csv"), row.names = FALSE)

      # 📊 **Volcano Plot**
      if (nrow(res) > 0) {
        p_volcano <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +  
          geom_point(aes(color = case_when(
            log2FoldChange < -high_log2fc_threshold & padj <= high_padj_threshold ~ "Low L2FC",
            log2FoldChange > high_log2fc_threshold & padj <= high_padj_threshold  ~ "High L2FC",
            TRUE ~ "Neutral"
          )), alpha = 0.7) +
          scale_color_manual(
            values = c("Low L2FC" = "#156082", "High L2FC" = "#990000", "Neutral" = "grey"),
            name = "",  # Removes legend title
            drop = FALSE  # Ensures all categories always appear
          ) +
          xlab("log2 Fold Change") +
          ylab("-log10 Adjusted p-value") +
          ggtitle(comparison_name) +
          theme_minimal() +
          theme(
            panel.grid = element_blank(),  # Removes gridlines
            plot.title = element_text(hjust = 0.5, face = "bold"),  # Centers and bolds title
            legend.position = "right",  # Keeps legend readable
            legend.text = element_text(size = 10)  # Makes legend text larger
          )

        # Save volcano plot
        ggsave(paste0("deseq/plots/volcano_", comparison_name, ".png"), plot = p_volcano, width = 8, height = 6, dpi = 300)

        # 📊 **MA Plot using DESeq2's built-in function**
        ma_plot_path <- paste0("deseq/plots/ma_", comparison_name, ".png")

        png(ma_plot_path, width = 8, height = 6, units = "in", res = 300)
        plotMA(res_deseq, main = paste("MA Plot:", comparison_name), ylim = c(-40, 40))
        dev.off()
      }
    }
  }
}

message("✅ All pairwise comparisons saved successfully.")
