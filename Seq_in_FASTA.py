from Bio import SeqIO

def count_sequences(fasta_file):
    count = 0
    with open(fasta_file, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            count += 1
    return count

if __name__ == "__main__":
    fasta_file = "MangoMeta_lncRNA/min_DELs.fasta"  # Replace with the path to your FASTA file
    num_sequences = count_sequences(fasta_file)
    print(f"Number of sequences in the file: {num_sequences}")
