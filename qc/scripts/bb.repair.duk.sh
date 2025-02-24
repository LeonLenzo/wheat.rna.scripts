#!/bin/bash

# Path to repair.sh script
repair_script="/usr/share/bbmap/repair.sh"

# Input directory containing raw reads
input_dir="./reads/raw"

# Loop through all R1 files in the input directory
for fastq in "$input_dir"/*-R1.fastq.gz; do
    sample_ID=$(basename "$fastq" -R1.fastq.gz)
    echo "Running repair.sh for sample: $sample_ID"

    # Repair command with increased memory allocation
    java -ea -Xmx60g -cp /usr/share/java/bbmap.jar jgi.SplitPairsAndSingles rp \
        in1="$input_dir/${sample_ID}-R1.fastq.gz" \
        in2="$input_dir/${sample_ID}-R2.fastq.gz" \
        out1="$input_dir/${sample_ID}-R1.repaired.fastq.gz" \
        out2="$input_dir/${sample_ID}-R2.repaired.fastq.gz" \
        outs="$input_dir/${sample_ID}.singletons.fastq.gz" \
        qin=33 tossbrokenreads nullifybrokenquality

    echo "Running bbduk for sample: $sample_ID"

    # BBDuk command
    java -ea -Xmx60g -cp /usr/share/java/bbmap.jar jgi.BBDuk \
        in1="$input_dir/${sample_ID}-R1.repaired.fastq.gz" \
        in2="$input_dir/${sample_ID}-R2.repaired.fastq.gz" \
        out1="$input_dir/${sample_ID}-R1.trimmed.fastq.gz" \
        out2="$input_dir/${sample_ID}-R2.trimmed.fastq.gz" \
        ref="./adapters.fa" \
        ktrimright=t \
        k=27 \
        hdist=1 \
        edist=0 \
        qtrim=rl \
        trimq=20 \
        minlength=40 \
        trimbyoverlap=t \
        minoverlap=24 \
        ordered=t \
        showspeed=t \
        qin=33 \
        stats="$input_dir/${sample_ID}.stats.txt"
done
