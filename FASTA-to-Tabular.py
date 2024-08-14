from Bio import SeqIO

def fasta_to_tabular(input_fasta, output_tsv):
    with open(output_tsv, 'w') as tsv_file:
        for record in SeqIO.parse(input_fasta, "fasta"):
            tsv_file.write(f"{record.id}\t{str(record.seq)}\n")

if __name__ == "__main__":
    input_fasta = "MangoMeta_lncRNA/mango_TA4_proteins.fasta"  # Replace with your input FASTA file path
    output_tsv = "output.tsv"    # Replace with your desired output TSV file path
    fasta_to_tabular(input_fasta, output_tsv)
