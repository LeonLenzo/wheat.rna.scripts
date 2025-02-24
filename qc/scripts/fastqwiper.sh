#!/bin/bash

# Input and output directories
input_dir="./3KO"
output_dir="./wiper"

# Docker image for FastqWiper
docker_image="mazzalab/fastqwiper"

# Number of cores to use for FastqWiper
num_cores=16

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

docker run --rm -ti 
--name fastqwiper 
-v "$input_dir YOUR_LOCAL_PATH_TO_DATA_FOLDER:/fastqwiper/data" mazzalab/fastqwiper paired 8 sample 33 ACGTN 500000

# Process all paired-end files in the input directory
for r1 in "$input_dir"/*-R1.fastq.gz; do
    # Identify the corresponding R2 file and sample name
    r2="${r1/-R1.fastq.gz/-R2.fastq.gz}"
    sample_name=$(basename "$r1" -R1.fastq.gz)

    # Ensure the paired file exists
    if [[ -f "$r2" ]]; then
        echo "Processing paired files: $sample_name"

        # Run FastqWiper Docker container for the paired files
        docker run --rm \
            -v "$input_dir:/fastqwiper/data" \
            -v "$output_dir:/fastqwiper/output" \
            "$docker_image" paired \
            "$num_cores" "$sample_name" 33 ACGTN 500000

        # Print the status
        if [[ $? -eq 0 ]]; then
            echo "Successfully processed: $sample_name"
        else
            echo "Error processing: $sample_name"
        fi
    else
        echo "Warning: No matching R2 file for $sample_name. Skipping."
    fi
done
