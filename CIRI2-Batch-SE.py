import os
import glob
import subprocess

# ==================== SET YOUR PATHS HERE ====================
INPUT_DIR  = "/run/media/joydeep/Expansion/Sorghum_circRNA/Drought/Cleaned_FASTQ/SE/"
OUTPUT_DIR = "/run/media/joydeep/Blue_Drive/SBI/Drought/CIRI2_Results_New/SE/"

REF_FASTA  = "/run/media/joydeep/Expansion/Sorghum_circRNA/BWA-MEM/sbi.fasta"
GTF_ANNO   = "/run/media/joydeep/Expansion/Sorghum_circRNA/GTF_Files/SbicircRNA_merged.gtf"

THREADS    = "10"
# ============================================================


def run_process():

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Change 1: Search for simple .fastq.gz files
    fastq_files = sorted(glob.glob(os.path.join(INPUT_DIR, "*.fastq.gz")))

    if not fastq_files:
        print(f"No FASTQ files found in {INPUT_DIR}")
        return

    print(f"[INFO] Found {len(fastq_files)} samples")

    for f1 in fastq_files:

        # Change 2: Get base name by stripping extension
        filename = os.path.basename(f1)
        base = filename.replace(".fastq.gz", "")

        print(f"\n>>> Processing {base}")

        # ---------- Sample-specific output directory ----------
        sample_dir = os.path.join(OUTPUT_DIR, base)
        os.makedirs(sample_dir, exist_ok=True)

        sam_file    = os.path.join(sample_dir, f"{base}.sam")
        bwa_log     = os.path.join(sample_dir, f"{base}_bwa.log")
        ciri_out    = os.path.join(sample_dir, f"{base}_ciri.txt")
        ciri_as_out = os.path.join(sample_dir, f"{base}_CIRI-AS.txt")

        # ---------- Step 1: BWA-MEM (Single End) ----------
        print("   [1/3] BWA-MEM alignment (Single-End)")
        bwa_cmd = (
            f"bwa mem -t {THREADS} -T 19 "
            f"{REF_FASTA} {f1} "
            f"1> {sam_file} 2> {bwa_log}"
        )
        subprocess.run(bwa_cmd, shell=True, check=True)

        # ---------- Step 2: CIRI2 ----------
        print("   [2/3] Running CIRI2")
        ciri_cmd = (
            f"CIRI.pl -I {sam_file} "
            f"-O {ciri_out} "
            f"-F {REF_FASTA} "
            f"-A {GTF_ANNO} "
            f"-T {THREADS}"
        )
        subprocess.run(ciri_cmd, shell=True, check=True)

        # ---------- Step 3: CIRI-AS ----------
        print("   [3/3] Running CIRI-AS")
        ciri_as_cmd = (
            f"CIRI-AS.pl "
            f"-S {sam_file} "
            f"-C {ciri_out} "
            f"-O {ciri_as_out} "
            f"-F {REF_FASTA} "
            f"-A {GTF_ANNO} "
            f"-D yes"
        )
        subprocess.run(ciri_as_cmd, shell=True, check=True)

        # ---------- Cleanup ----------
        if os.path.exists(sam_file):
            os.remove(sam_file)

    print("\n[DONE] CIRI2 + CIRI-AS pipeline completed successfully")


if __name__ == "__main__":
    run_process()