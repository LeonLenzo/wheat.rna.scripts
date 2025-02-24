import csv

# File paths
input_csv = "mapping/gene_to_uniprot_mapping.csv"  # Your original mapping file
output_csv = "mapping/gene_to_uniprot_expanded.csv"  # Expanded output file

# Expand the CSV so each UniProt ID gets its own row
def expand_csv(input_csv, output_csv):
    expanded_rows = []
    
    with open(input_csv, "r") as infile:
        reader = csv.reader(infile)
        header = next(reader)  # Read header
        
        for row in reader:
            gene_id, uniprot_ids = row[0], row[1]
            if uniprot_ids and uniprot_ids != "No UniProt IDs found":
                for uniprot_id in uniprot_ids.split(","):
                    expanded_rows.append([gene_id, uniprot_id.strip()])
            else:
                expanded_rows.append([gene_id, "No UniProt ID"])  # Keep unmapped entries
    
    with open(output_csv, "w", newline="") as outfile:
        writer = csv.writer(outfile)
        writer.writerow(["GeneID", "UniProt_ID"])  # Write header
        writer.writerows(expanded_rows)

# Run transformation
expand_csv(input_csv, output_csv)
print(f"Expanded Gene-to-UniProt mapping saved to '{output_csv}'")
