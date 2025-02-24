#!/bin/bash

# Input file containing gene IDs (one per line)
INPUT_FILE="/home/leonl/wheat.rna/mapping/IWGSC/v1.genes.txt"

# Output CSV file
OUTPUT_FILE="gene_mappings.csv"

# API base URL
ENSEMBL_API="https://rest.ensembl.org"

# Specify species (Ensuring wheat genes are mapped correctly)
SPECIES="triticum_aestivum"

# Check if input file exists
if [[ ! -f "$INPUT_FILE" ]]; then
    echo "Error: Input file '$INPUT_FILE' not found!"
    exit 1
fi

# Print header for CSV output
echo "Original_Gene_ID,Mapped_Gene_ID,Assembly" > "$OUTPUT_FILE"

# Read each gene ID from the input file and query the Ensembl REST API
while IFS= read -r GENE_ID; do
    # Query Ensembl API for gene mapping and assembly version
    RESPONSE=$(curl -s "${ENSEMBL_API}/lookup/id/${GENE_ID}?species=${SPECIES}&content-type=application/json")

    # Extract the mapped gene ID and assembly version using jq
    MAPPED_GENE_ID=$(echo "$RESPONSE" | jq -r '.id // "Not_Found"')
    ASSEMBLY_VERSION=$(echo "$RESPONSE" | jq -r '.assembly_name // "Unknown_Assembly"')

    # Write to output CSV
    echo "${GENE_ID},${MAPPED_GENE_ID},${ASSEMBLY_VERSION}" >> "$OUTPUT_FILE"

    # Optional: Show progress
    echo "Processed: $GENE_ID -> $MAPPED_GENE_ID (Assembly: $ASSEMBLY_VERSION)"

    # Sleep to avoid overloading the API (adjust as needed)
    sleep 0.2

done < "$INPUT_FILE"

echo "Gene mapping complete! Results saved to: $OUTPUT_FILE"
