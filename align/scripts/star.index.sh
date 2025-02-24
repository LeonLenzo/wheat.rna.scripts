#!/bin/bash

# Set paths
GENOME_FASTA="indices/GCF_018294505.1/GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_genomic.fna"  # Change to your genome FASTA file
GTF_FILE="indices/GCF_018294505.1/genomic.gff"    # Change to your GTF annotation file
OUTPUT_DIR="STAR_genome_index"
THREADS=20  # Adjust as needed

# Create output directory
mkdir -p ${OUTPUT_DIR}

STAR --runThreadN 20 \
     --runMode genomeGenerate \
     --genomeDir ${OUTPUT_DIR} \
     --genomeFastaFiles ${GENOME_FASTA} \
     --sjdbGTFfile ${GTF_FILE} \
     --sjdbOverhang 149 \
     --limitGenomeGenerateRAM 600000000 \ # 100GB
     --genomeChrBinNbits 16


echo "STAR genome index has been built!"
