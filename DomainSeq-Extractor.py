#!/usr/bin/env python3

from Bio import SeqIO

# ==== USER INPUT ====
# Path to your multi-FASTA file
fasta_file = "D:\\ST_Protein Motif\\STP.fasta"

# Path to your coordinate table (tab-separated)
coord_file = "D:\\ST_Protein Motif\\STP_CORD.txt"

# Output file for domain sequences
output_file = "D:\\ST_Protein Motif\\STP_Domain.fasta"

# ====================

# Load full sequences into a dictionary
seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

# Read coordinates
domain_coords = []
with open(coord_file, "r") as f:
    next(f)  # skip header
    for line in f:
        line = line.strip()
        if not line:
            continue
        protein_id, coords = line.split("\t")
        start, end = coords.split("-")
        domain_coords.append((protein_id, int(start), int(end)))

# Extract domains and write to new FASTA
with open(output_file, "w") as out:
    for protein_id, start, end in domain_coords:
        if protein_id not in seq_dict:
            print(f"Warning: {protein_id} not found in FASTA.")
            continue

        full_seq = seq_dict[protein_id].seq
        # Python is 0-based index, but coordinates are usually 1-based
        domain_seq = full_seq[start - 1 : end]
        # Write domain to file
        out.write(f">{protein_id}_{start}-{end}\n")
        out.write(f"{domain_seq}\n")

print(f"Domain sequences written to {output_file}")
