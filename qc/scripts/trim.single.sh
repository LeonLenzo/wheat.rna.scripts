#!/bin/bash

### Trim Raw Paired Reads ### 

# Input directory and adapter file
input_dir="./reads/raw"
adapters="./adapters.fa"

for fastq in "$input_dir"/*-R1.fastq.gz; do
    sample_ID=$(basename "$fastq" -R1.fastq.gz)

    bbduk.sh \
        in1="$input_dir/${sample_ID}-R1.fastq.gz" \
        in2="$input_dir/${sample_ID}-R2.fastq.gz" \
        out1="$input_dir/${sample_ID}-R1.trimmed.fastq.gz" \
        out2="$input_dir/${sample_ID}-R2.trimmed.fastq.gz" \
        ref="$adapters" \
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
