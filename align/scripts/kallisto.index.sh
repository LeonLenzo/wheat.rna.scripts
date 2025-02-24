#!/bin/bash

cDNA_file="indices/GCF_018294505.1/cds_from_genomic.fna"
index_file="indices/GCF_018294505.1/cds_from_genomic.idx"

kallisto index -i "$index_file" "$cDNA_file"