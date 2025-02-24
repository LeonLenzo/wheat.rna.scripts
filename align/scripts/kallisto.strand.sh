#!/bin/bash

# User-defined variables (modify these as needed)
KALLISTO_INDEX="indices/GCF_018294505.1/cds_from_genomic.idx"   # Path to kallisto index
OUT_DIR="reads/test/out"                                        # Output directory for results
SEQ_DIR="reads/test"                                            # Directory containing sequencing files
SAMPLE_NAME="inplanta-SN15-016"                                 # Sample name (prefix for files)

# Define input file paths
READ1="${SEQ_DIR}/${SAMPLE_NAME}-R1.trimmed.fastq.gz"
READ2="${SEQ_DIR}/${SAMPLE_NAME}-R2.trimmed.fastq.gz"

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"

# Check if input files exist
if [[ ! -f "$READ1" || ! -f "$READ2" ]]; then
    echo "Error: Input sequence files not found: $READ1 or $READ2"
    exit 1
fi

echo "Running kallisto quantification for different strandedness options..."

# Run kallisto quantification for different strandedness options
kallisto quant -i ${KALLISTO_INDEX} -o ${OUT_DIR}/${SAMPLE_NAME}.un ${READ1} ${READ2}
kallisto quant -i ${KALLISTO_INDEX} -o ${OUT_DIR}/${SAMPLE_NAME}.rf ${READ1} ${READ2} --rf-stranded
kallisto quant -i ${KALLISTO_INDEX} -o ${OUT_DIR}/${SAMPLE_NAME}.fr ${READ1} ${READ2} --fr-stranded

echo "Extracting abundance values and computing strandedness..."

# Extract abundance values and compute strandedness
paste ${OUT_DIR}/${SAMPLE_NAME}.fr/abundance.tsv \
      ${OUT_DIR}/${SAMPLE_NAME}.rf/abundance.tsv \
      ${OUT_DIR}/${SAMPLE_NAME}.un/abundance.tsv | \
cut -f1,4,9,14 | \
awk 'BEGIN{sum1=0;sum2=0;sum3=0}{sum1+=$2;sum2+=$3;sum3+=$4}END{print sum1, sum2, sum3}' > ${OUT_DIR}/${SAMPLE_NAME}.libtypetesting

echo "Determining library type..."

# Determine library type
less ${OUT_DIR}/${SAMPLE_NAME}.libtypetesting | \
awk '{print $2/$1, $3/$1, $3/$2}' | \
awk '{if($1<0.3 && $3>3) print "stranded"; else if($1>3 && $2>3) print "reverse"; else print "unstranded"}' \
> ${OUT_DIR}/${SAMPLE_NAME}.libtype

# Output result
echo "Library type determination complete. See ${OUT_DIR}/${SAMPLE_NAME}.libtype for results."
cat ${OUT_DIR}/${SAMPLE_NAME}.libtype
