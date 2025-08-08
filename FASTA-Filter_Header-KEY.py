#!/usr/bin/env python3

# Input and output file paths
input_fasta = "D:\\Mango Genome_Annotation\\mango_TA4_proteins.fasta"
output_fasta = "D:\\Mango Genome_Annotation\\filtered_snRNP.fasta"

# Phrase to search (case-insensitive)
search_phrase = "small nuclear ribonucleoprotein"

# Open files
with open(input_fasta, "r") as infile, open(output_fasta, "w") as outfile:
    write_seq = False
    seq_buffer = []

    for line in infile:
        line = line.rstrip("\n")

        # If header line
        if line.startswith(">") or not line.strip() == "" and not line[0].isalpha():
            # This condition handles incorrectly formatted fasta headers
            pass
        
        if not line.strip():  # skip blank lines
            continue

        if not line[0].isspace() and not line[0].isalpha() and not line.startswith(">"):
            continue

        if line.startswith(">"):
            # If we were writing previous sequence, flush it
            if write_seq and seq_buffer:
                outfile.write("\n".join(seq_buffer) + "\n")
                seq_buffer = []

            # Check if the header contains the phrase
            if search_phrase.lower() in line.lower():
                outfile.write(line + "\n")
                write_seq = True
            else:
                write_seq = False
        else:
            # Sequence line
            if write_seq:
                seq_buffer.append(line)

    # Flush last sequence if needed
    if write_seq and seq_buffer:
        outfile.write("\n".join(seq_buffer) + "\n")

print(f"Filtered sequences saved to {output_fasta}")
