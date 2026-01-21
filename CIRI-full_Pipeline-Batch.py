#!/usr/bin/env python3

import os
import subprocess
import re
import sys

# ---------------- USER SETTINGS ---------------- #

FASTQ_DIR = "/home/pgb-lab/Documents/Sorghum_circRNA/Drought/Cleaned_FASTQ"
REFERENCE = "/home/pgb-lab/Documents/Sorghum_circRNA/BWA-MEM/sbi.fasta"
ANNOTATION = "/home/pgb-lab/Documents/Sorghum_circRNA/Drought/DeNovo_Assembled/Sbi_merged_circ.gtf"
OUTDIR = "/home/pgb-lab/Documents/Sorghum_circRNA/Drought/CIRI-full_Output"
THREADS = 10

# ---------------------------------------------- #

# Create only the TOP-LEVEL output directory
os.makedirs(OUTDIR, exist_ok=True)

# Collect FASTQ files
fastq_files = [f for f in os.listdir(FASTQ_DIR) if f.endswith(".fastq.gz")]
if not fastq_files:
    sys.exit("ERROR: No fastq.gz files found")

# Pair reads
pairs = {}

for f in fastq_files:
    m = re.match(r"(.+?)_(1|2)\.fastq\.gz$", f)
    if not m:
        print(f"Skipping unrecognized file: {f}")
        continue

    sample, read = m.groups()
    pairs.setdefault(sample, {})[read] = f

# Run CIRI-full per sample
for sample, reads in sorted(pairs.items()):

    if not ("1" in reads and "2" in reads):
        print(f"Skipping incomplete pair for sample: {sample}")
        continue

    r1 = reads["1"]
    r2 = reads["2"]

    sample_outdir = os.path.join(OUTDIR, sample)

    # IMPORTANT: do NOT pre-create this directory
    if os.path.exists(sample_outdir):
        print(f"Skipping {sample}: output directory already exists")
        continue

    cmd = [
        "CIRI-full", "Pipeline",
        "-1", os.path.join(FASTQ_DIR, r1),
        "-2", os.path.join(FASTQ_DIR, r2),
        "-r", REFERENCE,
        "-d", sample_outdir,
        "-o", sample,
        "-t", str(THREADS),
        "-0"
    ]

    if ANNOTATION:
        cmd.extend(["-a", ANNOTATION])

    print(f"\n▶ Running CIRI-full for sample: {sample}")
    print("  ", " ".join(cmd))

    subprocess.run(cmd, check=True)

print("\n✅ All samples processed successfully")
