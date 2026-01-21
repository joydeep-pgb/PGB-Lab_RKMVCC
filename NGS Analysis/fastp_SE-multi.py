import subprocess
import logging
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm

# ==================== CONFIGURATION ====================
# Set your paths here
INPUT_DIR  = Path("/run/media/joydeep/One_HDD/Sorghum_MetaDEG/Drought/FASTQ/SE")
OUTPUT_DIR = Path("/run/media/joydeep/Blue_Drive/Circ_Trimmed/Drought/SE/")
REPORT_DIR = OUTPUT_DIR / "fastp-reports"
LOG_FILE   = "fastp_se_analysis.log"

# Performance Settings
NUM_WORKERS = 10      # Number of files to process in parallel
THREADS_PER_JOB = 2   # Threads fastp uses per file
# =======================================================

def setup_logging():
    logging.basicConfig(
        filename=LOG_FILE,
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def run_fastp_se(fastq_in):
    """Processes a single FASTQ file."""
    # Create sample name by removing extensions
    sample_name = fastq_in.name.split('.')[0]
    
    # Define outputs
    out_fastq = OUTPUT_DIR / f"{sample_name}_trimmed.fastq.gz"
    html_report = REPORT_DIR / f"{sample_name}_report.html"
    json_report = REPORT_DIR / f"{sample_name}_report.json"

    # Construct fastp command for Single-End
    cmd = [
        "fastp",
        "-i", str(fastq_in),
        "-o", str(out_fastq),
        "-h", str(html_report),
        "-j", str(json_report),
        "--thread", str(THREADS_PER_JOB),
        "--max_len1", "100",
        "--length_required", "100",
        "-R", f"Fastp SE report for {sample_name}"
    ]

    try:
        logging.info(f"Processing: {sample_name}")
        # Run process and capture output for logging if it fails
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        logging.info(f"Successfully processed: {sample_name}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error processing {sample_name}: {e.stderr}")

def main():
    # Ensure directories exist
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    
    setup_logging()

    # Find all gzipped fastq files
    fastq_files = list(INPUT_DIR.glob("*.fastq.gz"))
    
    if not fastq_files:
        print(f"Error: No .fastq.gz files found in {INPUT_DIR}")
        return

    print(f"Found {len(fastq_files)} files. Starting parallel Single-End trimming...")

    # Execute jobs in parallel with a progress bar
    with ProcessPoolExecutor(max_workers=NUM_WORKERS) as executor:
        list(tqdm(executor.map(run_fastp_se, fastq_files), total=len(fastq_files), desc="Trimming Progress"))

    print(f"\nProcessing complete. Detailed log: {LOG_FILE}")

if __name__ == "__main__":
    main()