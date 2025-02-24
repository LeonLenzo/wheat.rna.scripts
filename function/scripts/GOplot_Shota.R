# Load Required Libraries
library(ggplot2)
library(dplyr)
library(scales)
library(viridis)  # Install if needed: install.packages("viridis")

# Load GO Enrichment Results
go_results <- read.csv("GO/deseq2_Tween_vs_SN15_GO_enrichment.csv")  # Update filename if needed

# Clean GO Categories and Fix Whitespace Issues
go_results <- go_results %>%
  mutate(
    GO_Category = case_when(
      GO_Category == "MF" ~ "Molecular Function",
      GO_Category == "BP" ~ "Biological Process",
      GO_Category == "CC" ~ "Cellular Component",
      TRUE ~ GO_Category
    ),
    GO_Description = trimws(GO_Description)  # Remove leading/trailing spaces
  )

# Calculate Dynamic Axis Limits
x_min <- min(go_results$Sum_log2FoldChange, na.rm = TRUE)
x_max <- max(go_results$Sum_log2FoldChange, na.rm = TRUE)
x_breaks <- scales::pretty_breaks(n = 5)(c(x_min, x_max))  # Adjusted breaks

y_min <- min(go_results$p.adjust, na.rm = TRUE)
y_max <- max(go_results$p.adjust, na.rm = TRUE)
y_breaks <- scales::log_breaks()(c(y_min, y_max))  # Log-transformed breaks

# Generate Bubble Plot with **No Borders**
p <- ggplot(go_results, aes(
  x = Sum_log2FoldChange, 
  y = p.adjust, 
  fill = Count, 
  size = Count, 
  label = GO_Description)) +

  # **Use shape=16 (solid circles with NO borders)**
  geom_point(alpha = 0.9, shape = 16, aes(colour = GO_Category)) +  

  # Dynamic x-axis scaling
  scale_x_continuous(limits = c(x_min, x_max), breaks = x_breaks) +

  # Dynamic y-axis scaling (log10 transformation + inversion)
  scale_y_continuous(trans = scales::log10_trans(), 
                     limits = c(y_min, y_max),
                     breaks = y_breaks) +
  scale_y_reverse() +  # **Invert y-axis separately**

  # Styling and modifications
  scale_fill_viridis_c(option = "C", name = "Gene Count (%)") +
  scale_color_manual(values = c("#88CCEE","#117733","#332288"), 
                     labels = c("Biological Process", "Molecular Function", "Cellular Component")) +

  # Keep only the vertical dashed line at x = 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +

  labs(x = expression(Sigma * " log"[2] * " Fold Change"), 
       y = "Adjusted p-value",
       colour = "GO Category",
       size = "# Genes") +

  theme_minimal(base_size = 14) +  # **Use minimal theme to avoid hidden overrides**
  theme(
    axis.line = element_line(size = 1, colour = "black"),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.text = element_text(size=12),
    axis.title = element_text(size=14),
    legend.key = element_blank()  # Remove legend background
  )

# Save Bubble Plot
ggsave("GO_Enrichment_BubblePlot.pdf", plot = p, width = 10, height = 6)

# Print message indicating completion
print("✅ GO Enrichment Bubble Plot saved as 'GO_Enrichment_BubblePlot.pdf'")
