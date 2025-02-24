import time
import pandas as pd
from Bio import Entrez
import xml.etree.ElementTree as ET

# Set your email for NCBI Entrez
Entrez.email = "leon.lenzo.1@gmail.com"  # Replace with your actual email

# Load LOC numbers from CSV file
input_csv = "LOC.test.csv"  # Replace with your file
output_csv = "LOC_XM_XP_mappings.csv"
debug_output = "debug_full_reports.xml"  # Save full reports for debugging

df = pd.read_csv(input_csv)
loc_list = df.iloc[:, 0].tolist()  # Assuming LOC numbers are in the first column

# Function to query NCBI Gene for structured data using esummary
def fetch_gene_summary(loc_id):
    """Fetch structured summary XML for a given LOC Gene ID from NCBI."""
    gene_id = loc_id.replace("LOC", "")  # Remove LOC prefix
    try:
        handle = Entrez.esummary(db="gene", id=gene_id, retmode="xml")
        return handle.read().decode("utf-8")  # ✅ FIX: Decode bytes to string
    except Exception as e:
        return f"Error: {e}"

# Function to extract XM → XP mappings from XML
def extract_mappings(xml_text):
    """Parse XML summary and extract XM_* → XP_* mappings."""
    mappings = []
    
    try:
        root = ET.fromstring(xml_text)

        for docsum in root.findall(".//DocumentSummary"):
            # First, try extracting from OtherAliases
            for item in docsum.findall(".//OtherAliases"):
                aliases = item.text
                if aliases:
                    pairs = [line.split() for line in aliases.split(";") if "XM_" in line and "XP_" in line]
                    for pair in pairs:
                        if len(pair) >= 2:
                            mappings.append((pair[0], pair[1]))  # (XM, XP)

            # Alternative approach: Check if mappings exist in RNA/Protein fields
            for rna in docsum.findall(".//RNA"):
                xm_accession = rna.find("Accession")
                if xm_accession is not None and "XM_" in xm_accession.text:
                    for protein in docsum.findall(".//Protein"):
                        xp_accession = protein.find("Accession")
                        if xp_accession is not None and "XP_" in xp_accession.text:
                            mappings.append((xm_accession.text, xp_accession.text))

    except Exception as e:
        print(f"⚠️ XML Parsing Error: {e}")
    
    return mappings

# Process each LOC and extract mappings
results = []
for loc in loc_list:
    print(f"Processing {loc}...")
    xml_report = fetch_gene_summary(loc)

    if "Error" in xml_report:
        print(f"⚠️ Issue retrieving summary for {loc}: {xml_report}")
        continue  # Skip to the next LOC if the report is invalid

    xm_xp_mappings = extract_mappings(xml_report)

    if not xm_xp_mappings:
        print(f"⚠️ No XM → XP mappings found for {loc}.")
        with open(debug_output, "a") as debug_file:
            debug_file.write(f"\n----- DEBUG REPORT for {loc} -----\n")
            debug_file.write(xml_report)  # Save XML report for debugging
        continue  # Skip to next LOC if no mappings

    for xm, xp in xm_xp_mappings:
        print(f"Found Mapping: {xm} → {xp}")  # Log extracted mapping
        results.append([loc, xm, xp])

    time.sleep(1)  # Avoid NCBI rate limits

# Convert results to DataFrame and save
columns = ["LOC_ID", "XM_Accession", "XP_Accession"]
df_results = pd.DataFrame(results, columns=columns)

if df_results.empty:
    print("⚠️ No results were extracted. Check the debug XML outputs.")

df_results.to_csv(output_csv, index=False)
print(f"✅ LOC to XM/XP mappings saved to {output_csv}")
