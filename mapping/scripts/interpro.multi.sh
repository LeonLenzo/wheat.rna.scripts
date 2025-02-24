#!/bin/bash

INPUT_DIR="$HOME/interproscan-5.72-103.0/data/split_fasta"
OUTPUT_DIR="$HOME/interproscan-5.72-103.0/data/jobs/interproscan_results"
mkdir -p "$OUTPUT_DIR"

for file in "$INPUT_DIR"/*.fa; do
    # Ensure file exists (prevent issue when no files match the pattern)
    if [[ ! -f "$file" ]]; then
        echo "Warning: No FASTA files found in $INPUT_DIR!"
        exit 1
    fi

    output_file="$OUTPUT_DIR/$(basename "$file" .fa).tsv"
    echo "Processing $file -> $output_file"
    
docker run --rm \
    -v "$HOME/interproscan-5.72-103.0/data:/opt/interproscan/data" \
    -v "$HOME/interproscan-5.72-103.0/data/split_fasta:/data/split_fasta" \
    -v "$HOME/interproscan-5.72-103.0/data/jobs/interproscan_results:/jobs/interproscan_results" \
    interpro/interproscan:5.72-103.0 \
    --input "/data/split_fasta/$(basename "$file")" \
    --output-file "/jobs/interproscan_results/$(basename "$output_file")" \
    --cpu 16 \
    --goterms \
    --iprlookup \
    --formats TSV


done
