from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def tabular_to_fasta(input_tsv, output_fasta):
    records = []
    with open(input_tsv, 'r') as tsv_file:
        for line in tsv_file:
            line = line.strip()
            if line:
                parts = line.split('\t')
                if len(parts) != 2:
                    raise ValueError("Each line in the input TSV must have exactly two columns")
                identifier, sequence = parts
                record = SeqRecord(Seq(sequence), id=identifier, description="")
                records.append(record)
    
    with open(output_fasta, 'w') as fasta_file:
        SeqIO.write(records, fasta_file, "fasta")

if __name__ == "__main__":
    input_tsv = "output.tsv"  # Replace with your input TSV file path
    output_fasta = "MangoMeta_lncRNA/mango_TA4_proteins.fasta"  # Replace with your desired output FASTA file path
    tabular_to_fasta(input_tsv, output_fasta)
