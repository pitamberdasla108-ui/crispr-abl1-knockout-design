#!/usr/bin/env python3
"""fetch_abl1.py — Download the human ABL1 reference mRNA + CDS from NCBI."""
import os
from Bio import Entrez, SeqIO

Entrez.email = "pitamber.das.351@my.csun.edu"  # NCBI asks for an email
HERE = os.path.dirname(os.path.abspath(__file__))
os.makedirs(os.path.join(HERE,"data"), exist_ok=True)

# ABL1 RefSeq: NM_005157 is the canonical transcript (isoform a, 1b)
ACC = "NM_005157"

print(f"Fetching {ACC} (human ABL1) from NCBI...")
handle = Entrez.efetch(db="nucleotide", id=ACC, rettype="gb", retmode="text")
rec = SeqIO.read(handle, "genbank")
handle.close()

print(f"\nID: {rec.id}")
print(f"Description: {rec.description}")
print(f"Length: {len(rec.seq)} bp")

# find the CDS (coding sequence) and key features
cds = [f for f in rec.features if f.type == "CDS"]
if cds:
    c = cds[0]
    print(f"\nCDS location: {c.location}")
    print(f"CDS length: {len(c.location)} bp -> {len(c.location)//3} codons (~{len(c.location)//3-1} aa)")
    gene = c.qualifiers.get("gene", ["?"])[0]
    print(f"Gene: {gene}")

# save the full sequence + CDS
SeqIO.write(rec, os.path.join(HERE,"data","ABL1_NM_005157.gb"), "genbank")
with open(os.path.join(HERE,"data","ABL1_mRNA.fasta"),"w") as f:
    f.write(f">{rec.id} {rec.description}\n{str(rec.seq)}\n")
if cds:
    cds_seq = cds[0].location.extract(rec.seq)
    with open(os.path.join(HERE,"data","ABL1_CDS.fasta"),"w") as f:
        f.write(f">ABL1_CDS\n{str(cds_seq)}\n")
    print(f"\nFirst 90 bp of CDS (start of coding region):")
    print(f"  {str(cds_seq[:90])}")
    print(f"Starts with ATG (start codon)? {str(cds_seq[:3])=='ATG'}")

print("\nSaved: data/ABL1_NM_005157.gb, ABL1_mRNA.fasta, ABL1_CDS.fasta")
