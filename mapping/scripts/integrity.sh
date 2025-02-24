#!/bin/bash

# Define directories
input_dir="."                # Directory containing FASTQ files
log_dir="./seqkit_results"   # Directory containing integrity log files
output_dir="./cleaned_fastq" # Directory to save cleaned FASTQ files
removed_log_dir="./removed_logs" # Directory to save logs of removed reads

# Create necessary directories if they don't exist
mkdir -p "$output_dir" "$removed_log_dir"

# Process each FASTQ file
for fastq in "$input_dir"/*.fastq.gz; do
    # Extract file name and construct paths
    file_name=$(basename "$fastq") # Get only the file name
    integrity_log="$log_dir/${file_name}.integrity.log" # Path to the corresponding integrity log
    corrupt_reads_log="$removed_log_dir/${file_name%.fastq.gz}_corrupt_headers.txt"
    output_fastq="$output_dir/$file_name"

    # Debug: Show paths being used
    echo "Processing $fastq"
    echo "Looking for integrity log: $integrity_log"

    # Check if the integrity log exists
    if [[ ! -f "$integrity_log" ]]; then
        echo "Integrity log for $file_name not found. Skipping..."
        continue
    fi

    # Extract headers of corrupt reads
    echo "Extracting corrupt read headers..."
    grep "Error in record" "$integrity_log" | awk '{print $NF}' > "$corrupt_reads_log"

    # Check if any corrupt reads were identified
    if [[ ! -s "$corrupt_reads_log" ]]; then
        echo "No corrupt reads found in $file_name. Skipping file..."
        continue
    fi

    # Remove corrupt reads using seqkit
    echo "Removing corrupt reads from $file_name..."
    seqkit grep -v -f "$corrupt_reads_log" "$fastq" | gzip > "$output_fastq"

    # Log results
    echo "Cleaned file saved to $output_fastq"
    echo "Headers of removed reads saved to $corrupt_reads_log"
done

echo "Processing completed. Cleaned files are in $output_dir, and logs are in $removed_log_dir."
