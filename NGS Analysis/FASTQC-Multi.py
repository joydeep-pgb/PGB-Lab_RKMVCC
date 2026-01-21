import os
import subprocess
import logging
from pathlib import Path
from multiprocessing import Pool, cpu_count
from tqdm import tqdm  # Install via: pip install tqdm

# ================= CONFIGURATION =================
# Set your paths and parameters here
INPUT_FOLDER = "/run/media/joydeep/One_HDD/Sorghum_MetaDEG/Drought/FASTQ/SE/"
OUTPUT_FOLDER = "/run/media/joydeep/One_HDD/Sorghum_MetaDEG/Drought/FASTQ/Reports/"
NUM_WORKERS = 4  # Number of parallel processes
# =================================================

def setup_logging(output_dir):
    """Logs to a file in the output directory."""
    log_path = Path(output_dir) / 'fastqc_analysis.log'
    logging.basicConfig(
        filename=log_path,
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        filemode='w' # 'w' overwrites log each time, 'a' appends
    )

def run_fastqc(args_tuple):
    """The task performed by each worker."""
    fastq_path, output_dir = args_tuple
    
    # Building the command as a list is safer than a string
    command = ["fastqc", "-o", str(output_dir), str(fastq_path)]
    
    try:
        # We capture output to keep it from cluttering the tqdm progress bar
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        return True, fastq_path.name
    except subprocess.CalledProcessError as e:
        return False, f"{fastq_path.name}: {e.stderr}"

def main():
    # 1. Prepare Paths
    input_dir = Path(INPUT_FOLDER)
    output_dir = Path(OUTPUT_FOLDER)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 2. Setup Log
    setup_logging(output_dir)

    # 3. Find Files
    # This finds .fastq.gz, .fq.gz, .fastq, and .fq
    fastq_files = [f for f in input_dir.iterdir() if f.name.endswith(('.fastq.gz', '.fq.gz', '.fastq', '.fq'))]
    
    if not fastq_files:
        print(f"Error: No FASTQ files found in {INPUT_FOLDER}")
        return

    # 4. Process with Progress Bar
    tasks = [(f, output_dir) for f in fastq_files]
    
    print(f"Starting FastQC on {len(fastq_files)} files...")
    print(f"Using {NUM_WORKERS} parallel workers.")

    success_count = 0
    failure_count = 0

    with Pool(NUM_WORKERS) as pool:
        # imap_unordered allows the progress bar to update as soon as ANY file is done
        for success, message in tqdm(pool.imap_unordered(run_fastqc, tasks), total=len(tasks), desc="FastQC Progress"):
            if success:
                logging.info(f"Successfully processed: {message}")
                success_count += 1
            else:
                logging.error(f"Error processing {message}")
                failure_count += 1

    # 5. Final Summary
    print("\n--- Analysis Completed ---")
    print(f"Successfully processed: {success_count}")
    print(f"Failures: {failure_count}")
    print(f"Logs saved to: {output_dir / 'fastqc_analysis.log'}")

if __name__ == "__main__":
    main()