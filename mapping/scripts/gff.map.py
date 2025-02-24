import gzip

def extract_xm_xp_mapping(gff_file, output_file):
    xm_xp_map = {}

    # Open GFF file (handle both compressed and uncompressed formats)
    open_func = gzip.open if gff_file.endswith(".gz") else open
    with open_func(gff_file, "rt") as f:
        for line in f:
            if line.startswith("#"):  # Skip header lines
                continue

            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue  # Skip malformed lines

            feature_type = cols[2]  # e.g., "mRNA", "CDS"
            attributes = cols[8]  # Column 9 contains the attributes

            # Extract XM ID from mRNA features
            if feature_type == "mRNA":
                if "ID=" in attributes and "XM_" in attributes:
                    xm_id = attributes.split("ID=")[1].split(";")[0]
                    xm_xp_map[xm_id] = None  # Initialize with no XP yet

            # Extract XP ID from CDS features
            elif feature_type == "CDS":
                if "Parent=" in attributes and "protein_id=" in attributes:
                    parent_id = attributes.split("Parent=")[1].split(";")[0]
                    protein_id = attributes.split("protein_id=")[1].split(";")[0]
                    if parent_id in xm_xp_map:
                        xm_xp_map[parent_id] = protein_id

    # Write XM → XP mapping to output file
    with open(output_file, "w") as out_f:
        out_f.write("XM_Transcript_ID\tXP_Protein_ID\n")
        for xm, xp in xm_xp_map.items():
            out_f.write(f"{xm}\t{xp if xp else 'No XP found'}\n")

    print(f"Mapping completed. Output saved to: {output_file}")

# Example usage
gff_path = "indices/GCF_018294505.1/genomic.gff"  # Update with your GFF file path
output_path = "mapping/XF2XP.tsv"
extract_xm_xp_mapping(gff_path, output_path)
