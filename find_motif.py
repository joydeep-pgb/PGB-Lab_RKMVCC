#!/usr/bin/env python3

from Bio import SeqIO
import re

# Your PROSITE pattern as a regex:
# Translate:
# [LIVFYW]-x-[SAGTN]-E-[LIVFMTN]-x5-[RS]-[GRKYA]-x7-[LIVMFST]-x-[LIVWF]-x2-G-x-[LIVFM]-[LIVGTFY]
pattern = r"[LIVFYW].[SAGTN]E[LIVFMTN].{5}[RS][GRKYA].{7}[LIVMFST].[LIVWF]..G.[LIVFM][LIVGTFY]"

# Compile regex for speed
motif = re.compile(pattern)

# FASTA file
fasta_file = "G:\Bioinformatics_Data\mango_sugar_transproter\Final Protein Prop\Motif_Reanalysis\STP.fasta"

# Scan each sequence
for record in SeqIO.parse(fasta_file, "fasta"):
    seq = str(record.seq)
    for match in motif.finditer(seq):
        start = match.start() + 1  # 1-based indexing
        end = match.end()
        print(f"{record.id}\t{start}-{end}\t{match.group()}")
