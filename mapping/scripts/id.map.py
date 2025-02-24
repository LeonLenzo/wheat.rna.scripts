import os
import csv
import time
import threading
from Bio import Entrez
import xml.etree.ElementTree as ET

# Configure Entrez
Entrez.email = "leon.lenzo.1@gmail.com"
Entrez.tool = "GeneToUniProtMapper"

# Input and output file paths
input_file = "CS_CDS_GeneIDs.txt"  # Plain text file with Gene IDs (one per line)
output_file = "gene_to_uniprot_mapping.csv"

# Set the number of concurrent threads (adjust as needed)
max_threads = 5
thread_lock = threading.Lock()

# Function to download XML data for a given Gene ID
def download_xml(gene_id):
    try:
        handle = Entrez.efetch(db="gene", id=gene_id, rettype="xml")
        xml_data = handle.read()
        handle.close()
        return xml_data
    except Exception as e:
        print(f"Error downloading XML for Gene ID {gene_id}: {e}")
        return None

# Function to parse UniProt annotations from XML data
def parse_xml_and_get_uniprot(xml_data):
    try:
        root = ET.fromstring(xml_data)
        uniprot_ids = []

        for comment in root.findall(".//Gene-commentary"):
            heading = comment.find("Gene-commentary_heading")
            if heading is not None and heading.text == "UniProtKB":
                for source in comment.findall(".//Other-source"):
                    dbtag = source.find("Other-source_src/Dbtag")
                    if dbtag is not None and dbtag.find("Dbtag_db").text == "UniProtKB/TrEMBL":
                        uniprot_id = dbtag.find("Dbtag_tag/Object-id/Object-id_str")
                        if uniprot_id is not None:
                            uniprot_ids.append(uniprot_id.text)

        return ",".join(uniprot_ids) if uniprot_ids else "No UniProt IDs found"
    except Exception as e:
        print(f"Error parsing XML: {e}")
        return "Error"

# Function to process a single gene
def process_gene(gene_id, results):
    xml_data = download_xml(gene_id)
    if xml_data:
        uniprot_ids = parse_xml_and_get_uniprot(xml_data)

        # Print the result to stdout
        print(f"Gene ID: {gene_id} -> UniProt IDs: {uniprot_ids}")

        with thread_lock:
            results.append((gene_id, uniprot_ids))

    time.sleep(0.2)  # Throttle requests to avoid API limits

# Read GeneIDs from input text file
def read_gene_ids(file_path):
    with open(file_path, "r") as txtfile:
        return [line.strip() for line in txtfile if line.strip()]

# Write mapping results to output CSV
def write_results(results, file_path):
    with open(file_path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["GeneID", "UniProt_IDs"])
        writer.writerows(results)

# Main function
def main():
    gene_ids = read_gene_ids(input_file)
    results = []
    threads = []

    for gene_id in gene_ids:
        while threading.active_count() > max_threads:
            time.sleep(0.1)  # Wait for an available thread slot

        thread = threading.Thread(target=process_gene, args=(gene_id, results))
        thread.start()
        threads.append(thread)

    for thread in threads:
        thread.join()

    # Save results
    write_results(results, output_file)
    print(f"Gene-to-UniProt mapping saved to '{output_file}'")

if __name__ == "__main__":
    main()
