import pandas as pd

def table_to_fasta():
    table_file = "top_hits/top.csv"  # Hardcoded input file
    fasta_file = "top_hits/top.faa"  # Hardcoded output file
    delimiter = ","  # Hardcoded delimiter (change to "\t" for TSV)
    
    df = pd.read_csv(table_file, delimiter=delimiter)
    
    if "Title" not in df.columns or "Sequence" not in df.columns:
        raise ValueError("Table must contain 'Title' and 'Sequence' columns.")
    
    with open(fasta_file, "w") as fasta_out:
        for _, row in df.iterrows():
            fasta_out.write(f">{row['Title']}\n{row['Sequence']}\n")
    
    print(f"FASTA file saved to {fasta_file}")

if __name__ == "__main__":
    table_to_fasta()
