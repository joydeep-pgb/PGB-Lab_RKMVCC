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

    Parameters:
    file_path (str): Path to the file containing SRR accession IDs.

    Returns:
    list: List of SRR accession IDs.
    """
    try:
        with open(file_path, 'r') as file:
            return [line.strip() for line in file if line.strip()]
    except FileNotFoundError:
        logging.error(f"File not found: {file_path}")
        return []

def prefetch_sra_file(srr: str, output_dir: str, max_retries: int = 3) -> bool:
    """
    Download a single SRA file using prefetch.

    Parameters:
    srr (str): SRR accession ID.
    output_dir (str): Directory to store the downloaded file.
    max_retries (int): Maximum number of retry attempts.

    Returns:
    bool: True if the download was successful, False otherwise.
    """
    command = ["prefetch", "--output-directory", output_dir, srr]
    logging.info(f"Running command: {' '.join(command)}")
    
    for attempt in range(max_retries):
        try:
            result = subprocess.run(command, check=True, capture_output=True, text=True)
            logging.info(f"Output: {result.stdout}")  # Log output
            logging.info(f"Successfully downloaded {srr}")
            return True
        except subprocess.CalledProcessError as e:
            logging.error(f"Attempt {attempt + 1} failed for {srr}: {e}")
            logging.error(f"Error Output: {e.stderr}")  # Log error output
            if attempt < max_retries - 1:
                time.sleep(2)  # Wait before retrying
            else:
                logging.error(f"Max retries reached for {srr}. Skipping.")
                return False

def fastq_dump_single_file(srr: str, sra_dir: str, fastq_dir: str) -> bool:
    """
    Convert a single SRA file to FASTQ files for paired-end reads using fastq-dump.

    Parameters:
    srr (str): SRR accession ID.
    sra_dir (str): Directory where the downloaded SRA files are stored.
    fastq_dir (str): Directory to store the FASTQ files.

    Returns:
    bool: True if the conversion was successful, False otherwise.
    """
    sra_file = os.path.join(sra_dir, srr, f"{srr}.sra")  # Point to the correct .sra file
    if not os.path.isfile(sra_file):
        logging.warning(f"SRA file for {srr} not found in {sra_dir}. Skipping.")
        return False
    
    command = ["fastq-dump", "--split-files", "--gzip", "--outdir", fastq_dir, sra_file]
    logging.info(f"Running command: {' '.join(command)}")
    try:
        subprocess.run(command, check=True)
        logging.info(f"Successfully converted {srr} to FASTQ files")
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to convert {srr} to FASTQ files: {e}")
        return False

def prefetch_sra_files_parallel(srr_accessions: List[str], output_dir: str, max_workers: int = 4) -> None:
    """
    Download SRA files in parallel using ThreadPoolExecutor.

    Parameters:
    srr_accessions (list): List of SRR accession IDs.
    output_dir (str): Directory to store the downloaded files.
    max_workers (int): Maximum number of threads to use.

    Returns:
    None
    """
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(prefetch_sra_file, srr, output_dir): srr for srr in srr_accessions}
        for future in as_completed(futures):
            srr = futures[future]
            try:
                success = future.result()
                if success:
                    logging.info(f"Completed download for {srr}")
                else:
                    logging.error(f"Failed to download {srr}")
            except Exception as e:
                logging.error(f"Unexpected error for {srr}: {e}")

def fastq_dump_paired_end_parallel(srr_accessions: List[str], sra_dir: str, fastq_dir: str, max_workers: int = 4) -> None:
    """
    Convert SRA files to FASTQ files in parallel using ThreadPoolExecutor.

    Parameters:
    srr_accessions (list): List of SRR accession IDs.
    sra_dir (str): Directory where the downloaded SRA files are stored.
    fastq_dir (str): Directory to store the FASTQ files.
    max_workers (int): Maximum number of threads to use.

    Returns:
    None
    """
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(fastq_dump_single_file, srr, sra_dir, fastq_dir): srr for srr in srr_accessions}
        for future in as_completed(futures):
            srr = futures[future]
            try:
                success = future.result()
                if success:
                    logging.info(f"Completed conversion for {srr}")
                else:
                    logging.error(f"Failed to convert {srr}")
            except Exception as e:
                logging.error(f"Unexpected error for {srr}: {e}")

def main():
    # Path to the file containing SRR accession IDs
    srr_file_path = "/mnt/d/LAB_Data/Phaseolus_DEG/Phaseolus_Root/SUS/SRA_List.txt"
    
    # Directory to download SRA files to
    sra_dir = "/mnt/d/LAB_Data/Phaseolus_DEG/Phaseolus_Root/SUS/SRA/"
    
    # Directory to store FASTQ files
    fastq_dir = "/mnt/d/LAB_Data/Phaseolus_DEG/Phaseolus_Root/SUS/FASTQ/Drought_T75/"

    # Create directories if they don't exist
    os.makedirs(sra_dir, exist_ok=True)
    os.makedirs(fastq_dir, exist_ok=True)

    # Read SRR accessions from file
    srr_accessions = read_srr_accessions(srr_file_path)
    if not srr_accessions:
        logging.error("No SRR accessions found. Exiting.")
        return

    # Download SRA files in parallel
    logging.info("Starting parallel download of SRA files...")
    prefetch_sra_files_parallel(srr_accessions, sra_dir, max_workers=4)

    # Convert SRA files to FASTQ files in parallel
    logging.info("Starting parallel conversion of SRA files to FASTQ...")
    fastq_dump_paired_end_parallel(srr_accessions, sra_dir, fastq_dir, max_workers=4)

if __name__ == "__main__":
    main()
