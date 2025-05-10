import subprocess
import os
import logging
import time
from typing import List
from concurrent.futures import ThreadPoolExecutor, as_completed

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def read_srr_accessions(file_path: str) -> List[str]:
    """
    Read SRR accessions from a file.
    """
    try:
        with open(file_path, 'r') as file:
            return [line.strip() for line in file if line.strip()]
    except FileNotFoundError:
        logging.error(f"File not found: {file_path}")
        return []

def prefetch_sra_file(srr: str, output_dir: str, max_retries: int = 3) -> bool:
    """
    Download a single SRA file using prefetch with retries.
    """
    command = ["prefetch", "--output-directory", output_dir, srr]
    logging.info(f"Running command: {' '.join(command)}")
    
    for attempt in range(max_retries):
        try:
            result = subprocess.run(command, check=True, capture_output=True, text=True)
            logging.info(f"Output: {result.stdout}")
            return True
        except subprocess.CalledProcessError as e:
            logging.error(f"Attempt {attempt + 1} failed for {srr}: {e}")
            logging.error(f"Error Output: {e.stderr}")
            if attempt < max_retries - 1:
                time.sleep(2)
            else:
                logging.error(f"Max retries reached for {srr}. Skipping.")
                return False

def fastq_dump_single_file(srr: str, sra_dir: str, fastq_dir: str) -> bool:
    """
    Convert a single SRA file to FASTQ using fastq-dump.
    """
    sra_file = os.path.join(sra_dir, srr, f"{srr}.sra")
    if not os.path.isfile(sra_file):
        logging.warning(f"SRA file for {srr} not found. Skipping.")
        return False
    
    command = ["fastq-dump", "--split-files", "--gzip", "--outdir", fastq_dir, sra_file]
    logging.info(f"Running command: {' '.join(command)}")
    try:
        subprocess.run(command, check=True)
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to convert {srr}: {e}")
        return False

def process_srr(srr: str, sra_dir: str, fastq_dir: str, max_retries: int = 3) -> bool:
    """
    Handle download and conversion for a single SRR accession.
    """
    # Download SRA file
    download_success = prefetch_sra_file(srr, sra_dir, max_retries)
    if not download_success:
        return False
    
    # Convert to FASTQ
    conversion_success = fastq_dump_single_file(srr, sra_dir, fastq_dir)
    return conversion_success

def main():
    # Configuration paths
    srr_file_path = "/home/pgb-lab/Documents/Pvul_SRA/Drought_Root.txt"
    sra_dir = "/home/pgb-lab/Documents/Pvul_SRA/SRA/"
    fastq_dir = "/home/pgb-lab/Documents/Pvul_SRA/FASTQ/"

    # Create directories if needed
    os.makedirs(sra_dir, exist_ok=True)
    os.makedirs(fastq_dir, exist_ok=True)

    # Read SRR accessions
    srr_accessions = read_srr_accessions(srr_file_path)
    if not srr_accessions:
        logging.error("No SRR accessions found. Exiting.")
        return

    # Process all SRRs in parallel with immediate conversion after download
    with ThreadPoolExecutor(max_workers=4) as executor:
        futures = {executor.submit(process_srr, srr, sra_dir, fastq_dir): srr for srr in srr_accessions}
        for future in as_completed(futures):
            srr = futures[future]
            try:
                success = future.result()
                if success:
                    logging.info(f"Successfully processed {srr}")
                else:
                    logging.error(f"Failed to process {srr}")
            except Exception as e:
                logging.error(f"Unexpected error processing {srr}: {e}")

if __name__ == "__main__":
    main()