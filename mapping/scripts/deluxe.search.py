import os
import csv
from Bio import Entrez
import xml.etree.ElementTree as ET

# Configure Entrez
Entrez.email = "leon.lenzo.1@gmail.com"
Entrez.tool = "XMtoXPMapper"

# Input and output file paths
input_file = "top_hits/gene_list.txt"  # Plain text file with XM IDs (one per line)
output_file = "top_hits/xm_to_xp_mapping.csv"
expanded_output_file = "top_hits/xm_to_xp_expanded.csv"  # Expanded table

# Function to download XML data for a given XM accession
def download_xml(xm_accession, output_path):
    try:
        print(f"Downloading XML for XM Accession: {xm_accession}")
        handle = Entrez.efetch(db="nuccore", id=xm_accession, rettype="xml")
        with open(output_path, "wb") as file:
            file.write(handle.read())
        handle.close()
        return True
    except Exception as e:
        print(f"Error downloading XML for XM Accession {xm_accession}: {e}")
        return False

# Function to parse XP annotations, gene description, and protein sequence from the XML file
def parse_xml_and_get_xp(xml_file_path):
    try:
        with open(xml_file_path, "rb") as file:
            tree = ET.parse(file)
            root = tree.getroot()

            xp_accessions = []
            gene_description = "Unknown"
            protein_sequence = "Not found"

            # Extract gene description
            description_element = root.find(".//Entrezgene_gene/Gene-ref/Gene-ref_desc")
            if description_element is not None:
                gene_description = description_element.text

            # Extract XP accessions and protein sequence
            for product in root.findall(".//Gene-commentary"):
                heading = product.find("Gene-commentary_heading")
                if heading is not None and heading.text == "RefSeq Proteins":
                    for sub_product in product.findall(".//Gene-commentary"):
                        # Get XP accession
                        accession = sub_product.find("Gene-commentary_accession")
                        if accession is not None and accession.text.startswith("XP_"):
                            xp_accessions.append(accession.text)
                        
                        # Get protein sequence
                        seq_element = sub_product.find(".//Seq-loc_whole/Seq-id/Seq-id_gi")
                        if seq_element is not None:
                            protein_sequence = seq_element.text  # Placeholder for actual sequence retrieval

            return xp_accessions, gene_description, protein_sequence
    except Exception as e:
        print(f"Error parsing XML: {e}")
        return [], "Error", "Error"

# Read XM accessions from input text file
def read_xm_accessions(file_path):
    with open(file_path, "r") as txtfile:
        return [line.strip() for line in txtfile if line.strip()]

# Write mapping results to output CSV
def write_results(results, file_path):
    with open(file_path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["XM_Accession", "XP_Accession(s)", "Gene_Description", "Protein_Sequence"])
        writer.writerows(results)

# Expand the CSV so each XP ID gets its own row
def expand_csv(input_csv, output_csv):
    expanded_rows = []
    with open(input_csv, "r") as infile:
        reader = csv.reader(infile)
        header = next(reader)
        for row in reader:
            xm_id, xp_ids, gene_desc, protein_seq = row
            for xp_id in xp_ids.split(","):
                expanded_rows.append([xm_id, xp_id.strip(), gene_desc, protein_seq])
    
    with open(output_csv, "w", newline="") as outfile:
        writer = csv.writer(outfile)
        writer.writerow(["XM_Accession", "XP_Accession", "Gene_Description", "Protein_Sequence"])
        writer.writerows(expanded_rows)

# Main function
def main():
    xm_accessions = read_xm_accessions(input_file)
    results = []

    for xm_accession in xm_accessions:
        xml_file_path = f"{xm_accession}.xml"
        if download_xml(xm_accession, xml_file_path):  # Download XML
            xp_accessions, gene_description, protein_sequence = parse_xml_and_get_xp(xml_file_path)  # Parse XP IDs
            xp_accession_str = ",".join(xp_accessions) if xp_accessions else "No XP IDs found"
            print(f"XM Accession {xm_accession}: XP Accessions -> {xp_accession_str}, Gene -> {gene_description}")  # Sanity check output
            results.append((xm_accession, xp_accession_str, gene_description, protein_sequence))
            os.remove(xml_file_path)  # Remove temporary file
        else:
            results.append((xm_accession, "Error", "Error", "Error"))

    # Save results
    write_results(results, output_file)
    print(f"XM-to-XP mapping saved to '{output_file}'")

    # Expand the results into separate rows
    expand_csv(output_file, expanded_output_file)
    print(f"Expanded XM-to-XP mapping saved to '{expanded_output_file}'")

if __name__ == "__main__":
    main()
