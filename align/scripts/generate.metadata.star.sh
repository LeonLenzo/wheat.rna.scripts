#!/bin/bash

# Directory containing STAR ReadsPerGene.out.tab files
INPUT_DIR="/home/leonl/wheat.rna/align/star/"
OUTPUT_FILE="/home/leonl/wheat.rna/align/star/metadata.csv"

# Get list of all ReadsPerGene.out.tab files
FILES=($(ls ${INPUT_DIR}/*ReadsPerGene.out.tab 2>/dev/null))

# Check if files exist
if [[ ${#FILES[@]} -eq 0 ]]; then
    echo "Error: No ReadsPerGene.out.tab files found in '$INPUT_DIR'."
    exit 1
fi

# Print header for metadata file
echo "sample,condition,path" > "$OUTPUT_FILE"

# Process each file
for file in "${FILES[@]}"; do
    # Extract sample name
    SAMPLE=$(basename "$file" | sed 's/_ReadsPerGene.out.tab//')

    # Assign condition based on sample name
    if [[ "$SAMPLE" =~ Tween ]]; then
        CONDITION="Tween"
    elif [[ "$SAMPLE" =~ SN15 ]]; then
        CONDITION="SN15"
    elif [[ "$SAMPLE" =~ Pf2 ]]; then
        CONDITION="Pf2"
    elif [[ "$SAMPLE" =~ 3KO ]]; then
        CONDITION="3KO"
    else
        CONDITION="unknown"  # Assign "unknown" if no match
    fi

    # Write metadata row
    echo "$SAMPLE,$CONDITION,$file" >> "$OUTPUT_FILE"
done

echo "Metadata file saved to: $OUTPUT_FILE"
