from Bio import SeqIO

def clean_fasta(input_file, output_file):
    ambiguous_characters = set("NRYWSMKHBVD")
    total_sequences = 0
    cleaned_sequences = 0
    total_removed_chars = 0

    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for seq_record in SeqIO.parse(infile, "fasta"):
            total_sequences += 1

            # Original sequence
            original_seq = str(seq_record.seq).upper()

            # Remove ambiguous characters
            cleaned_seq = ''.join([base for base in original_seq if base not in ambiguous_characters])
            removed_chars = len(original_seq) - len(cleaned_seq)
            total_removed_chars += removed_chars

            # Update sequence
            seq_record.seq = seq_record.seq.__class__(cleaned_seq)
            SeqIO.write(seq_record, outfile, "fasta")
            cleaned_sequences += 1

    print(f"Total sequences processed: {total_sequences}")
    print(f"Total sequences written: {cleaned_sequences}")
    print(f"Total ambiguous characters removed: {total_removed_chars}")

# Usage
input_file = "/media/pgb-lab/One_HDD/Mango_AS_rMATS/IsoformSwitch/IsoformNT.fasta"
output_file = "/media/pgb-lab/One_HDD/Mango_AS_rMATS/IsoformSwitch/IsoformNT_Cleaned.fasta"
clean_fasta(input_file, output_file)
