from Bio import SeqIO
import pandas as pd
import argparse

def fasta_to_table(fasta_file, output_file, output_format="csv"):
    records = []
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        title_modified = record.description + "$"  # Append $ at the end of the title
        records.append([title_modified, len(record.seq), str(record.seq)])
    
    df = pd.DataFrame(records, columns=["Title", "Length", "Sequence"])
    
    if output_format == "csv":
        df.to_csv(output_file, index=False)
    elif output_format == "tsv":
        df.to_csv(output_file, index=False, sep="\t")
    else:
        raise ValueError("Unsupported output format. Choose 'csv' or 'tsv'.")
    
    print(f"Table saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert FASTA to a tabular format.")
    parser.add_argument("fasta_file", help="Input FASTA file")
    parser.add_argument("output_file", help="Output table file (CSV or TSV)")
    parser.add_argument("--format", choices=["csv", "tsv"], default="csv", help="Output format (default: csv)")
    
    args = parser.parse_args()
    fasta_to_table(args.fasta_file, args.output_file, args.format)