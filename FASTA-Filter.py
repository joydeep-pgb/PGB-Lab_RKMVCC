from Bio import SeqIO

def filter_fasta(input_file, output_file):
    ambiguous_characters = set("NRYWSMKHBVD")
    total_sequences = 0
    filtered_sequences = 0

    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        sequences = SeqIO.parse(infile, "fasta")
        
        for seq_record in sequences:
            total_sequences += 1
            if not any(char in ambiguous_characters for char in seq_record.seq.upper()):
                SeqIO.write(seq_record, outfile, "fasta")
                filtered_sequences += 1

    print(f"Total sequences before filtering: {total_sequences}")
    print(f"Total sequences after filtering: {filtered_sequences}")

# Usage
input_file = "MangoMeta_lncRNA/FC_1/trans/GFvsRF.DEGs.fasta"
output_file = "MangoMeta_lncRNA/FC_1/trans/GFvsRF.DEGs.filtered_output.fasta"
filter_fasta(input_file, output_file)
