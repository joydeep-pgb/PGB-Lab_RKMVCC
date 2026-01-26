#!/usr/bin/env python3

import os
import subprocess
import re
import sys
import shutil

# ---------------- USER SETTINGS ---------------- #
FASTQ_DIR = "/run/media/joydeep/Expansion/Sorghum_circRNA/Drought/Cleaned_FASTQ/"
REFERENCE = "/run/media/joydeep/Expansion/Sorghum_circRNA/BWA-MEM/sbi.fasta"
ANNOTATION = "/run/media/joydeep/Expansion/Sorghum_circRNA/GTF_Files/SbicircRNA_merged.gtf"
OUTDIR = "/home/joydeep/Documents/SBI/Drought/CIRI-full_Output"
THREADS = 10
# ---------------------------------------------- #

os.makedirs(OUTDIR, exist_ok=True)

fastq_files = [f for f in os.listdir(FASTQ_DIR) if f.endswith(".fastq.gz")]
if not fastq_files:
    sys.exit("ERROR: No fastq.gz files found")

pairs = {}
for f in fastq_files:
    m = re.match(r"(.+?)_(1|2)\.fastq\.gz$", f)
    if not m:
        print(f"Skipping unrecognized file: {f}")
        continue
    sample, read = m.groups()
    pairs.setdefault(sample, {})[read] = f

for sample, reads in sorted(pairs.items()):
    if not ("1" in reads and "2" in reads):
        print(f"Skipping incomplete pair for sample: {sample}")
        continue

    r1 = reads["1"]
    r2 = reads["2"]
    sample_outdir = os.path.join(OUTDIR, sample)

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

    try:
        subprocess.run(cmd, check=True)
        
        # --- TARGETED CLEANUP ---
        # Define paths based on your provided structure
        files_to_remove = [
            os.path.join(sample_outdir, "sam", f"{sample}_ciri.sam"),
            os.path.join(sample_outdir, "sam", f"{sample}_ro.sam"),
            os.path.join(sample_outdir, "CIRI-full_output", f"{sample}_ro2.sam")
        ]
        
        for file_path in files_to_remove:
            if os.path.exists(file_path):
                print(f"🧹 Removing: {file_path}")
                os.remove(file_path)
            else:
                print(f"ℹ️ Note: {file_path} not found, skipping.")

        # Optional: Remove the 'sam' directory if it's now empty
        sam_dir = os.path.join(sample_outdir, "sam")
        if os.path.exists(sam_dir) and not os.listdir(sam_dir):
            os.rmdir(sam_dir)

    except subprocess.CalledProcessError as e:
        print(f"❌ Error processing sample {sample}: {e}")
        continue

print("\n✅ All samples processed and space cleared.")