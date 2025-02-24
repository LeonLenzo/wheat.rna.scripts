import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Directories containing the CSV files
filtered_directory = "deseq/high_l2fc_filtered_results"
unfiltered_directory = "deseq/raw_results"

# List all CSV files in both directories
filtered_files = [f for f in os.listdir(filtered_directory) if f.endswith(".csv")]
unfiltered_files = {f.replace(".csv", ""): os.path.join(unfiltered_directory, f) for f in os.listdir(unfiltered_directory) if f.endswith(".csv")}

# Dictionary to store data from each dataset
gene_list = set()

# Extract gene names from filtered datasets
for file in filtered_files:
    file_path = os.path.join(filtered_directory, file)
    df = pd.read_csv(file_path)
    
    # Standardize column names
    df.columns = df.columns.str.strip()
    
    if "gene" in df.columns:
        gene_list.update(df["gene"].tolist())
    else:
        print(f"Warning: 'gene' column not found in {file}")

gene_list = sorted(gene_list)  # Ensure consistent ordering

# Create a dataframe with gene names
merged_df = pd.DataFrame({"gene": gene_list})

# Extract log2FoldChange values from unfiltered datasets
for dataset_name, file_path in unfiltered_files.items():
    df = pd.read_csv(file_path)
    df.columns = df.columns.str.strip()  # Fix column spacing issues

    if "gene" in df.columns and "log2FoldChange" in df.columns:
        df = df.loc[:, ["gene", "log2FoldChange"]]  # Explicitly select columns by name
        df.rename(columns={"log2FoldChange": dataset_name}, inplace=True)
        
        # Merge with the main dataframe
        merged_df = pd.merge(merged_df, df, on="gene", how="left")
    else:
        print(f"Skipping {file_path}: Missing required columns.")

# Fill missing values with 'NA' for clarity
merged_df.fillna("NA", inplace=True)

# Save merged results
output_file = os.path.join(unfiltered_directory, "merged_gene_data.csv")
merged_df.to_csv(output_file, index=False)
print(f"Merged gene data saved to {output_file}")

# Load the merged gene data for heatmap
df = pd.read_csv(output_file)

# Convert "NA" to NaN and ensure all log2FoldChange values are numeric
df.replace("NA", float("nan"), inplace=True)
df.set_index("gene", inplace=True)
df = df.astype(float)  # Convert all dataset columns to float

# Transpose the DataFrame to switch axes
df = df.T

# Check if DataFrame is empty before plotting
if df.empty:
    print("Warning: No data available for heatmap.")
else:
    plt.figure(figsize=(20, 10))  # Increase figure size for better readability
    heatmap = sns.heatmap(df, cmap="coolwarm", center=0, annot=False, linewidths=0.5)
    
    # Customize heatmap appearance
    plt.title("Gene Expression Heatmap (log2FoldChange)")
    plt.xlabel("Genes")
    plt.ylabel("Datasets")
    plt.xticks(rotation=90, ha="right", fontsize=8)  # Reduce font size for clarity
    plt.yticks(fontsize=10)
    
    # Ensure colorbar is correctly assigned
    cbar = heatmap.collections[0].colorbar
    cbar.set_label("log2FoldChange")
    
    # Save the heatmap image
    heatmap_output_path = os.path.join(unfiltered_directory, "gene_expression_heatmap_landscape.png")
    plt.savefig(heatmap_output_path, dpi=300, bbox_inches='tight')
    print(f"Heatmap saved to {heatmap_output_path}")
    
    # Show plot
    plt.show()
