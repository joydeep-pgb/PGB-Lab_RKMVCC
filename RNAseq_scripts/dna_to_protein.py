import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def convert_to_protein(dna_sequence):
    dna_seq = Seq(dna_sequence)
    protein_seq = dna_seq.translate()
    return protein_seq

if len(sys.argv) != 3:
    print("Usage: python dna_to_protein.py input.fasta output.fasta")
    sys.exit(1)

input_fasta = sys.argv[1]
output_fasta = sys.argv[2]

protein_records = []

for record in SeqIO.parse(input_fasta, "fasta"):
    protein_sequence = convert_to_protein(str(record.seq))
    protein_record = SeqRecord(protein_sequence, id=record.id, description=record.description)
    protein_records.append(protein_record)

SeqIO.write(protein_records, output_fasta, "fasta")

print(f"Conversion complete. Protein sequences saved in {output_fasta}")
