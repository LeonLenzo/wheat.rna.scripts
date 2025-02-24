#!/bin/bash

# Script to run InterProScan using Docker

# Define variables
DOCKER_IMAGE="interpro/interproscan:5.72-103.0"  # Use a specific version
INPUT_FILE="./Triticum_aestivum_mace.PGSBv2.1.cdna.all.fa"  # Input file
OUTPUT_FILE="annotation_output.tsv"  # Output file
WORK_DIR=$(pwd)  # Current working directory
DATA_DIR="$WORK_DIR/interproscan-5.72-103.0/data"  # InterProScan data directory
LICENSED_DIR="$WORK_DIR/licensed"  # Directory for licensed components (optional)

# Check if the input file exists
if [[ ! -f "$WORK_DIR/$INPUT_FILE" ]]; then
    echo "Error: Input file '$INPUT_FILE' not found in $WORK_DIR."
    exit 1
fi

# Check if the data directory exists
if [[ ! -d "$DATA_DIR" ]]; then
    echo "Error: Data directory '$DATA_DIR' not found."
    exit 1
fi

# Run the Docker container
echo "Running InterProScan with Docker..."
docker run --rm \
    -v "$WORK_DIR:/data" \
    -v "$DATA_DIR:/opt/interproscan/data" \
    -v "$LICENSED_DIR:/opt/interproscan/licensed" \
    $DOCKER_IMAGE \
    interproscan.sh -i "/data/$INPUT_FILE" -f tsv -o "/data/$OUTPUT_FILE"

# Check if the output file was generated
if [[ -f "$WORK_DIR/$OUTPUT_FILE" ]]; then
    echo "InterProScan completed successfully."
    echo "Output file: $WORK_DIR/$OUTPUT_FILE"
else
    echo "Error: Output file not generated."
    exit 1
fi
