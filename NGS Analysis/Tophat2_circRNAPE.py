#!/usr/bin/env python

import os
import subprocess
import shutil
import sys

# --------------------------------------------------
# Paths (Adjust gtf_path and index_path as needed)
# --------------------------------------------------
input_dir = "/home/pgb-lab/Documents/Sorghum_circRNA/Drought/Cleaned_FASTQ"
output_dir = "/home/pgb-lab/Documents/Sorghum_circRNA/Drought/circRNA_Mapped/TopHat_Output/"
index_path = "/home/pgb-lab/Documents/Sorghum_circRNA/sbi_bowtie_Index/sbi_genome"
gtf_path = "/home/pgb-lab/Documents/Sorghum_circRNA/Sbicolor.gtf" 
summary_dir = os.path.join(output_dir, "Summary")

if not os.path.exists(output_dir): os.makedirs(output_dir)
if not os.path.exists(summary_dir): os.makedirs(summary_dir)

# --------------------------------------------------
# Detect FASTQ pairs
# --------------------------------------------------
fastq_files_1 = sorted([f for f in os.listdir(input_dir) if f.endswith("_1.fastq.gz")])

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

    sample_out_dir = os.path.join(output_dir, sample)
    stderr_log = os.path.join(summary_dir, sample + "_tophat.stderr.log")

    print("\n" + "="*70)
    print("RUNNING TOPHAT2 FOR SAMPLE: {}".format(sample))
    print("="*70)

    # TopHat2 Command (Optimized for CIRCexplorer2)
    tophat_cmd = [
        "tophat2",
        "-p", "10",
        "-G", gtf_path,
        "-o", sample_out_dir,
        "--fusion-search",
        "--keep-fasta-order",
        "--bowtie1",
        "--no-coverage-search",
        index_path,
        fq1, fq2
    ]

    # --------------------------------------------------
    # Execute with Real-Time Output to Terminal and File
    # --------------------------------------------------
    # stderr=subprocess.STDOUT merges error logs and progress logs into one stream
    proc = subprocess.Popen(
        tophat_cmd, 
        stdout=subprocess.PIPE, 
        stderr=subprocess.STDOUT, 
        bufsize=1, 
        universal_newlines=True
    )

    with open(stderr_log, "w") as log_file:
        # This loop reads the TopHat2 output line-by-line as it happens
        for line in proc.stdout:
            sys.stdout.write(line) # Prints to terminal
            log_file.write(line)   # Writes to the .log file
            sys.stdout.flush()     # Forces immediate display
            log_file.flush()

    proc.wait()

    if proc.returncode != 0:
        print("\n[ERROR] TopHat2 failed for sample: {}".format(sample))
        continue

    # --------------------------------------------------
    # Post-Processing
    # --------------------------------------------------
    final_bam = os.path.join(sample_out_dir, "accepted_hits.bam")
    tophat_log = os.path.join(sample_out_dir, "align_summary.txt")
    summary_file = os.path.join(summary_dir, sample + "_tophat_summary.txt")

    if os.path.exists(final_bam):
        print("\n[INFO] Indexing BAM for {}...".format(sample))
        subprocess.call(["samtools", "index", final_bam])
        
        if os.path.exists(tophat_log):
            shutil.copy(tophat_log, summary_file)
        
        print("[OK] Successfully finished: {}".format(sample))
    else:
        print("[ERROR] Mapping finished but accepted_hits.bam not found.")

print("\n[COMPLETE] All samples have been processed.")