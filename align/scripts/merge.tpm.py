import os
import pandas as pd
import glob

# Define the path to the parent directory containing Kallisto output directories
kallisto_dir = "quant/chinese_spring"  # Change this to your actual path

# Find all abundance.tsv files in subdirectories
abundance_files = glob.glob(os.path.join(kallisto_dir, "*", "abundance.tsv"))

# Dictionary to store TPM data
tpm_dict = {}

# Process each abundance.tsv file
for file in abundance_files:
    sample_name = os.path.basename(os.path.dirname(file))  # Extract sample name from directory
    df = pd.read_csv(file, sep='\t', usecols=["target_id", "tpm"])  # Load only required columns
    df.rename(columns={"tpm": sample_name}, inplace=True)  # Rename TPM column to sample name
    
    # Merge with existing data
    if not tpm_dict:
        tpm_dict["target_id"] = df["target_id"]
    
    tpm_dict[sample_name] = df[sample_name]

# Convert dictionary to DataFrame
tpm_merged_df = pd.DataFrame(tpm_dict)

# Save merged TPM table
output_file = os.path.join(kallisto_dir, "merged_tpm_counts.tsv")
tpm_merged_df.to_csv(output_file, sep='\t', index=False)

print(f"TPM counts merged and saved to {output_file}")
