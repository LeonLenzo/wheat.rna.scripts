import pandas as pd
import re

# Load your TSV file into a DataFrame
input_file = "rna.tsv"  # Replace with your actual file
df = pd.read_csv(input_file, sep="\t", dtype=str, low_memory=False)  # Fix dtype warning

# Define the column containing LOC numbers mixed with text
column_with_loc = "Description"  # Change this to match your file

# Function to extract LOC numbers using regex
def extract_loc(text):
    match = re.search(r"LOC\d+", str(text))  # Find LOC followed by numbers
    return match.group(0) if match else "Not Found"

# Apply the function to extract LOC numbers
df["LOC_ID"] = df[column_with_loc].apply(extract_loc)

# Save the cleaned data
output_file = "cleaned_with_loc.tsv"
df.to_csv(output_file, sep="\t", index=False)

print(f"Extracted LOC numbers saved to {output_file}")
