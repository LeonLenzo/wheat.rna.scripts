#!/bin/bash

# Define the directory containing DESeq2 results
SOURCE_DIR="/home/leonl/wheat.rna/deseq/kallisto/filtered_results"
DEST_DIR="/home/leonl/wheat.rna/deseq/combined/"
ALIGNMENT_METHOD="NCBI"  # Change this to the correct method name

# Ensure the destination directory exists
mkdir -p "$DEST_DIR"

# Loop through all files that start with "deseq"
for file in "$SOURCE_DIR"/deseq2*; do
    if [[ -f "$file" ]]; then
        # Get the new filename by replacing "deseq" with the alignment method
        new_file="$DEST_DIR/$(basename "$file" | sed "s/deseq/$ALIGNMENT_METHOD/")"

        # Rename (move) the file
        mv "$file" "$new_file"
        echo "Renamed: $file → $new_file"
    fi
done

echo "✅ All DESeq2 files renamed and moved successfully!"

