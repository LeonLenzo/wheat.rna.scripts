import os
from Bio import Entrez
from time import sleep

# Configure Entrez
Entrez.email = "leon.lenzo.1@gmail.com"
Entrez.tool = "GeneToUniProtMapper"

# Input and output file paths
input_file = "indices/GCF_018294505.1/gene_ids.txt"  # Plain text file with Gene IDs (one per line)
output_dir = "gene_xml_files"  # Directory to save XML files

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Function to fetch and save XML file for a Gene ID
def fetch_and_save_gene_xml(gene_id):
    try:
        print(f"Downloading XML for Gene ID: {gene_id}")
        handle = Entrez.efetch(db="gene", id=gene_id, rettype="xml")
        xml_data = handle.read().decode('utf-8')  # Decode bytes to string
        handle.close()
        
        # Save the XML file
        xml_file_path = os.path.join(output_dir, f"{gene_id}.xml")
        with open(xml_file_path, "w", encoding="utf-8") as xml_file:
            xml_file.write(xml_data)
        print(f"Saved XML for Gene ID {gene_id} to {xml_file_path}")
    except Exception as e:
        print(f"Error downloading XML for Gene ID {gene_id}: {e}")

# Read Gene IDs from input file
def read_gene_ids(file_path):
    with open(file_path, "r") as txtfile:
        return [line.strip() for line in txtfile if line.strip()]

# Main function to download XML files
def main():
    gene_ids = read_gene_ids(input_file)
    for gene_id in gene_ids:
        fetch_and_save_gene_xml(gene_id)
        sleep(0.01)  # Add a delay to prevent hitting API rate limits

if __name__ == "__main__":
    main()
