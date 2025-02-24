#!/bin/bash

# Directory containing STAR ReadsPerGene.out.tab files
INPUT_DIR="/home/leonl/wheat.rna/align/star/"
OUTPUT_FILE="/home/leonl/wheat.rna/align/star/counts_matrix.csv"

# Get list of all ReadsPerGene.out.tab files
FILES=($(ls ${INPUT_DIR}/*ReadsPerGene.out.tab 2>/dev/null))

# Check if files exist
if [[ ${#FILES[@]} -eq 0 ]]; then
    echo "Error: No ReadsPerGene.out.tab files found in '$INPUT_DIR'."
    exit 1
fi

# Extract sample names from filenames
SAMPLE_NAMES=()
for file in "${FILES[@]}"; do
    SAMPLE_NAMES+=("$(basename "$file" | sed 's/_ReadsPerGene.out.tab//')")
done

# Extract gene IDs from the first file (Column 1) and remove STAR headers (first 4 lines)
awk 'NR>4 {print $1}' "${FILES[0]}" > gene_ids.txt

# Extract Column 4 (Reverse Stranded Counts) from all files and store in temp files
for file in "${FILES[@]}"; do
    awk 'NR>4 {print $4}' "$file" > "$(basename "$file").counts"
done

# Use `paste` to merge gene IDs with count columns
paste -d ',' gene_ids.txt $(ls *.counts) > temp_counts.txt

# Create the final output with proper column headers
echo -n "GeneID" > "$OUTPUT_FILE"
for sample in "${SAMPLE_NAMES[@]}"; do
    echo -n ",$sample" >> "$OUTPUT_FILE"
done
echo "" >> "$OUTPUT_FILE"

# Append the count data
cat temp_counts.txt >> "$OUTPUT_FILE"

# Cleanup
rm temp_counts.txt gene_ids.txt *.counts

echo "Counts matrix saved to: $OUTPUT_FILE"
