#!/usr/bin/env python

import os
import subprocess
import shutil
import sys

# --------------------------------------------------
# Paths
# --------------------------------------------------
input_dir = "/home/pgb-lab/Documents/Sorghum_circRNA/Drought/Cleaned_FASTQ"
output_dir = "/home/pgb-lab/Documents/Sorghum_circRNA/Drought/circRNA_Mapped/STAR_Output/"
index_path = "/home/pgb-lab/Documents/Sorghum_circRNA/sbi_STAR_Index"
summary_dir = os.path.join(output_dir, "Summary")

# --------------------------------------------------
# Create directories
# --------------------------------------------------
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

if not os.path.exists(summary_dir):
    os.makedirs(summary_dir)

# --------------------------------------------------
# Detect FASTQ pairs
# --------------------------------------------------
fastq_files_1 = sorted(
    f for f in os.listdir(input_dir) if f.endswith("_1.fastq.gz")
)

if not fastq_files_1:
    sys.exit("ERROR: No *_1.fastq.gz files found")

# --------------------------------------------------
# Process samples
# --------------------------------------------------
for fastq1 in fastq_files_1:

    sample = fastq1.replace("_1.fastq.gz", "")
    fastq2 = fastq1.replace("_1.fastq.gz", "_2.fastq.gz")

    fq1 = os.path.join(input_dir, fastq1)
    fq2 = os.path.join(input_dir, fastq2)

    if not os.path.exists(fq2):
        print("[SKIP] Missing mate file for {}".format(sample))
        continue

    prefix = os.path.join(output_dir, sample + "_")
    final_bam = os.path.join(output_dir, sample + "_sorted.bam")
    summary_file = os.path.join(summary_dir, sample + "_STAR_summary.txt")
    stderr_log = os.path.join(summary_dir, sample + "_STAR.stderr.log")

    print("[INFO] Processing sample: {}".format(sample))

    # --------------------------------------------------
    # STAR command (circRNA-optimized)
    # --------------------------------------------------
    star_cmd = [
        "STAR",
        "--runThreadN", "10",
        "--genomeDir", index_path,
        "--readFilesIn", fq1, fq2,
        "--readFilesCommand", "zcat",
        "--outSAMtype", "BAM", "SortedByCoordinate",
        "--chimSegmentMin", "10",
        "--outFileNamePrefix", prefix
    ]

    # --------------------------------------------------
    # Run STAR (Python 2 compatible)
    # --------------------------------------------------
    err = open(stderr_log, "w")
    ret = subprocess.call(
        star_cmd,
        stdout=open(os.devnull, "w"),
        stderr=err
    )
    err.close()

    if ret != 0:
        print("[ERROR] STAR failed for {}".format(sample))
        print("        See log: {}".format(stderr_log))
        continue

    # --------------------------------------------------
    # Validate outputs
    # --------------------------------------------------
    star_bam = prefix + "Aligned.sortedByCoord.out.bam"
    star_log = prefix + "Log.final.out"

    if not os.path.exists(star_bam):
        print("[ERROR] BAM missing for {}".format(sample))
        continue

    os.rename(star_bam, final_bam)

    if os.path.exists(star_log):
        shutil.move(star_log, summary_file)
    else:
        print("[WARN] Log.final.out missing for {}".format(sample))

    # --------------------------------------------------
    # Index BAM
    # --------------------------------------------------
    ret = subprocess.call(["samtools", "index", final_bam])
    if ret != 0:
        print("[ERROR] samtools index failed for {}".format(sample))
        continue

    print("[OK] Completed {}".format(sample))

print("[DONE] STAR mapping, sorting, and indexing completed.")
