import subprocess
import os
import logging
import time

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def read_srr_accessions(file_path):
    """
    Read SRR accessions from a file.

    Parameters:
    file_path (str): Path to the file containing SRR accession IDs.

    Returns:
    list: List of SRR accession IDs.
    """
    with open(file_path, 'r') as file:
        return [line.strip() for line in file if line.strip()]

def prefetch_sra_files(srr_accessions, output_dir):
    """
    Download SRA files for a list of SRR accessions using prefetch.

    Parameters:
    srr_accessions (list): List of SRR accession IDs.
    output_dir (str): Directory to store the downloaded files.

    Returns:
    None
    """
    for srr in srr_accessions:
        command = ["prefetch", "--output-directory", output_dir, srr]
        logging.info(f"Running command: {' '.join(command)}")
        for attempt in range(3):  # Retry up to 3 times
            try:
                result = subprocess.run(command, check=True, capture_output=True, text=True)
                logging.info(f"Output: {result.stdout}")  # Log output
                logging.info(f"Successfully downloaded {srr}")
                break
            except subprocess.CalledProcessError as e:
                logging.error(f"Attempt {attempt + 1} failed for {srr}: {e}")
                logging.error(f"Error Output: {e.stderr}")  # Log error output
                time.sleep(2)  # Wait before retrying

def fastq_dump_paired_end(srr_accessions, sra_dir, fastq_dir):
    """
    Convert SRA files to FASTQ files for paired-end reads using fastq-dump.

    Parameters:
    srr_accessions (list): List of SRR accession IDs.
    sra_dir (str): Directory where the downloaded SRA files are stored.
    fastq_dir (str): Directory to store the FASTQ files.

    Returns:
    None
    """
    for srr in srr_accessions:
        sra_file = os.path.join(sra_dir, srr, f"{srr}.sra")  # Point to the correct .sra file
        if not os.path.isfile(sra_file):
            logging.warning(f"SRA file for {srr} not found in {sra_dir}. Skipping.")
            continue
        
        command = ["fastq-dump", "--split-files", "--gzip", "--outdir", fastq_dir, sra_file]
        logging.info(f"Running command: {' '.join(command)}")
        try:
            subprocess.run(command, check=True)
            logging.info(f"Successfully converted {srr} to FASTQ files")
        except subprocess.CalledProcessError as e:
            logging.error(f"Failed to convert {srr} to FASTQ files: {e}")

if __name__ == "__main__":
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

    # Download SRA files
    prefetch_sra_files(srr_accessions, sra_dir)

    # Convert SRA files to FASTQ files
    fastq_dump_paired_end(srr_accessions, sra_dir, fastq_dir)
