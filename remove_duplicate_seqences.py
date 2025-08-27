#!/usr/bin/env python3
"""
Remove duplicate sequences from a FASTA file based on sequence content, not headers.
Keeps the first occurrence of each unique sequence.
"""

def remove_duplicates(input_fasta, output_fasta):
    seen = set()
    with open(input_fasta, "r") as infile, open(output_fasta, "w") as outfile:
        seq_id = None
        seq = []

        for line in infile:
            line = line.strip()
            if line.startswith(">"):  # header line
                if seq_id and seq:  # write previous seq if unique
                    sequence = "".join(seq)
                    if sequence not in seen:
                        seen.add(sequence)
                        outfile.write(f"{seq_id}\n{sequence}\n")
                seq_id = line
                seq = []
            else:
                seq.append(line)

        # write the last sequence
        if seq_id and seq:
            sequence = "".join(seq)
            if sequence not in seen:
                seen.add(sequence)
                outfile.write(f"{seq_id}\n{sequence}\n")

if __name__ == "__main__":
    # 🔹 Set your input and output FASTA files here
    input_fasta = "D:\\Phaseolus_Auxin\\Phaseolus miRNA DB\\pvu-miRNA.txt"
    output_fasta = "D:\\Phaseolus_Auxin\\Phaseolus miRNA DB\\output_unique.fasta"

    remove_duplicates(input_fasta, output_fasta)
    print(f"Finished. Unique sequences written to {output_fasta}")
