from Bio import SeqIO # type: ignore
import os

# Configurations
input_fasta = "mapping/protein.faa"  # Your big FASTA file
output_dir = "mapping/split_fasta"  # Directory for smaller FASTA files
sequences_per_file = 1000  # Adjust based on memory availability

# Create output directory if not exists
os.makedirs(output_dir, exist_ok=True)

# Read and split the FASTA file
def split_fasta(input_fasta, output_dir, sequences_per_file):
    records = list(SeqIO.parse(input_fasta, "fasta"))
    total_files = (len(records) // sequences_per_file) + 1

    for i in range(total_files):
        chunk = records[i * sequences_per_file : (i + 1) * sequences_per_file]
        if not chunk:
            break
        output_file = os.path.join(output_dir, f"split_{i+1}.fa")
        SeqIO.write(chunk, output_file, "fasta")
        print(f"Created: {output_file}")

split_fasta(input_fasta, output_dir, sequences_per_file)
