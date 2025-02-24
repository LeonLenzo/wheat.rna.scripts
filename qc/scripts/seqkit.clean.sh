#!/bin/bash

# Input and output directories
input_dir="."
output_dir="./cleaned_fastq"
log_dir="./removed_logs"
temp_dir="./temp_chunks"
chunk_size="3G"  # Size of each chunk (adjust for memory constraints)

# Create output and log directories if they don't exist
mkdir -p "$output_dir" "$log_dir" "$temp_dir"

# Process each FASTQ file in the input directory
for fastq in "$input_dir"/*.fastq.gz; do
    # Extract the file name
    file_name=$(basename "$fastq")
    echo "Processing $file_name..."

    # Split the file into smaller chunks
    echo "Splitting $file_name into chunks of $chunk_size..."
    zcat "$fastq" | split -C "$chunk_size" - "$temp_dir/${file_name%.fastq.gz}_chunk_"

    # Process each chunk in parallel
    chunk_files=("$temp_dir/${file_name%.fastq.gz}_chunk_"*)
    echo "Processing ${#chunk_files[@]} chunks in parallel..."
    parallel --gnu -j 10 seqkit seq -v --remove-corrupt {} -o {}.cleaned ::: "${chunk_files[@]}"

    # Combine cleaned chunks
    echo "Combining cleaned chunks for $file_name..."
    cat "$temp_dir/${file_name%.fastq.gz}_chunk_"*.cleaned | gzip > "$output_dir/$file_name"

    # Extract and save removed reads
    echo "Logging removed reads for $file_name..."
    seqkit grep -j 10 -f <(seqkit seq -v --remove-corrupt "$fastq" -o /dev/null | grep "corrupt read" | awk '{print $NF}') \
        -v "$fastq" | gzip > "$log_dir/${file_name%.fastq.gz}_removed_reads.fastq.gz"

    # Cleanup temporary files for this FASTQ
    rm -f "$temp_dir/${file_name%.fastq.gz}_chunk_"*
    echo "Processing for $file_name completed. Cleaned file and logs saved."
done

# Cleanup temporary directory
rm -rf "$temp_dir"

echo "All files processed. Cleaned files are in $output_dir and removed reads are saved in $log_dir."
