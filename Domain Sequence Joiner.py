from Bio import SeqIO
from collections import defaultdict

# Input and output files
input_fasta = "G:\\Bioinformatics_Data\\mango_sugar_transproter\\Final Protein Prop\\Motif_Reanalysis\\Domain_Seq\\SWEET_Domain.fasta"
output_fasta = "G:\\Bioinformatics_Data\\mango_sugar_transproter\\Final Protein Prop\\Motif_Reanalysis\\Domain_Seq\\SWEET_Domain_merged.fasta"

# Dictionary to store merged sequences
seq_dict = defaultdict(str)

# Parse the FASTA and merge sequences
for record in SeqIO.parse(input_fasta, "fasta"):
    seq_dict[record.id] += str(record.seq)

# Write merged sequences to a new FASTA
with open(output_fasta, "w") as out_handle:
    for seq_id, sequence in seq_dict.items():
        out_handle.write(f">{seq_id}\n")
        # Wrap sequence lines at 60 characters (optional)
        for i in range(0, len(sequence), 60):
            out_handle.write(sequence[i:i+60] + "\n")

print(f"Merged FASTA written to: {output_fasta}")
