#!/bin/bash

# Define the inputs and output directories
idx_file="indices/GCF_018294505.1/cds_from_genomic.idx"
input_dir="/mnt/d/wheat.rna"
output_dir="./quant/cs_cds"
mkdir -p "$output_dir"

# Loop through all R1 read files in the input directory
for fastq in "$input_dir"/*-R1.fastq.gz; do

    # Define the sample ID
    sample_ID=$(basename "$fastq" -R1.fastq.gz)

    # Define the read files
    R1="$input_dir/${sample_ID}-R1.fastq.gz"
    R2="$input_dir/${sample_ID}-R2.fastq.gz"

    # Check if R1 and R2 exist
    if [[ ! -f "$R1" || ! -f "$R2" ]]; then
        echo "Error: Missing files for $sample_ID. Skipping..."
        continue
    fi

    # Define the output directory for the sample
    sample_output_dir="$output_dir/$sample_ID"
    mkdir -p "$sample_output_dir"

    # Run kallisto quant
    echo "Running kallisto quant for $sample_ID..."
    kallisto quant \
        -i "$idx_file" \
        -o "$sample_output_dir" \
        -b 100 \
        -t 18 \
        "$R1" "$R2"
    echo "Finished $sample_ID"

done

echo "All samples processed. Results are in $output_dir."

