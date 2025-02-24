from Bio import Entrez

# Set email for NCBI
Entrez.email = "leon.lenzo.1@gmail.com"

# Fetch a single Gene ID
gene_id = "34688808"
handle = Entrez.efetch(db="gene", id=gene_id, rettype="xml")
records = Entrez.read(handle)
handle.close()

print(records)
