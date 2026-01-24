import subprocess
import logging
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm

# ==================== CONFIGURATION ====================
# Define your paths here
INPUT_DIR  = Path("/run/media/joydeep/One_HDD/Sorghum_MetaDEG/Salt/FASTQ/")
OUTPUT_DIR = Path("/run/media/joydeep/Blue_Drive/SBI/Salt/Circ_Trimmed/")
REPORT_DIR = OUTPUT_DIR / "fastp-reports"
LOG_FILE   = "fastp_analysis.log"

# Performance Settings
NUM_WORKERS = 10     # Number of samples to process at once
THREADS_PER_JOB = 1  # Threads fastp uses per sample (e.g., 2 threads * 10 workers = 20 total CPUs used)
# =======================================================

def setup_logging():
    logging.basicConfig(
        filename=LOG_FILE,
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def run_fastp(fastq_r1):
    """Processes a single pair of FASTQ files."""
    # Derive filenames
    sample_name = fastq_r1.name.replace("_1.fastq.gz", "")
    fastq_r2 = fastq_r1.parent / fastq_r1.name.replace("_1.fastq.gz", "_2.fastq.gz")
    
    out_r1 = OUTPUT_DIR / f"{sample_name}_1.fastq.gz"
    out_r2 = OUTPUT_DIR / f"{sample_name}_2.fastq.gz"
    html_report = REPORT_DIR / f"{sample_name}_report.html"
    json_report = REPORT_DIR / f"{sample_name}_report.json"

    # Safety check: Ensure R2 exists
    if not fastq_r2.exists():
        logging.error(f"SKIPPED: Missing R2 file for {sample_name}")
        return

    # Construct command
    cmd = [
        "fastp",
        "-i", str(fastq_r1),
        "-I", str(fastq_r2),
        "-o", str(out_r1),
        "-O", str(out_r2),
        "-h", str(html_report),
        "-j", str(json_report),
        "--thread", str(THREADS_PER_JOB),
        "--max_len1", "100",
        "--max_len2", "100",
        "--length_required", "100",
        "-R", f"Fastp report for {sample_name}"
    ]

    try:
        logging.info(f"Processing: {sample_name}")
        # check=True will raise an error if fastp fails
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        logging.info(f"Successfully processed: {sample_name}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error processing {sample_name}: {e.stderr}")

def main():
    # Create directories if they don't exist
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    
    setup_logging()

    # Find all R1 files
    fastq_files_r1 = list(INPUT_DIR.glob("*_1.fastq.gz"))
    
    if not fastq_files_r1:
        print(f"Error: No files found in {INPUT_DIR}")
        return

    print(f"Found {len(fastq_files_r1)} samples. Starting parallel trimming...")

    # Execute in parallel
    with ProcessPoolExecutor(max_workers=NUM_WORKERS) as executor:
        list(tqdm(executor.map(run_fastp, fastq_files_r1), total=len(fastq_files_r1), desc="Trimming Progress"))

    print(f"\nTask Complete. Log saved to {LOG_FILE}")

if __name__ == "__main__":
    main()