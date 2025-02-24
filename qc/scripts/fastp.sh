#!/bin/bash

### Trim Raw Paired Reads ### 

# Get the input directory and adpater files from the command line argument
input_dir="."
adapters="adapters.fa"

for fastq in $input_dir/*-R1.fastq.gz; do

    # Extract sample IDs without extensions
    sample_ID=$(basename -s -R1.fastq.gz "$fastq")
    
    fastp \
    -i "$input_dir/${sample_ID}-R1.fastq.gz" \
    -I "$input_dir/${sample_ID}-R2.fastq.gz" \
    -o "$input_dir/${sample_ID}-R1.trimmed.fastq.gz" \
    -O "$input_dir/${sample_ID}-R2.trimmed.fastq.gz" \
    --thread 20 \
    -5 -3 -M 30 \
    --html ${sample_ID}.fastp.html --json ${sample_ID}.fastp.json \
    --verbose 

done
