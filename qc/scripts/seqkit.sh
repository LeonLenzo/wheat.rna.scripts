#!/bin/bash

# Directory containing FASTQ files
input_dir="."
output_dir="./seqkit_results"

# Create output directory for seqkit results if it doesn't exist
mkdir -p "$output_dir"

# Loop through each FASTQ file in the directory
for fastq in "$input_dir"/*.fastq.gz; do
    # Extract file name without path
    file_name=$(basename "$fastq")
    echo "Analyzing $file_name..."
    
    # Check for sequence and quality consistency
    integrity_log="${output_dir}/${file_name}.integrity.log"
    zcat "$fastq" | awk '{
        if (NR % 4 == 2) seq_length = length($0)
        if (NR % 4 == 0) {
            if (length($0) != seq_length) {
                print "Error in record:", NR / 4
                exit 1
            }
        }
    }' > "$integrity_log" 2>&1

    # Check if the integrity log contains errors
    if [[ -s "$integrity_log" ]]; then
        echo "Integrity issue detected in $file_name. Check log: $integrity_log"
        continue
    else
        echo "No issues found in $file_name."
    fi

    # Generate basic statistics for the file
    stats_file="${output_dir}/${file_name}.stats.txt"
    seqkit stats "$fastq" > "$stats_file"
    echo "Statistics saved to $stats_file"
done

echo "Seqkit analysis completed. Results are in $output_dir."
