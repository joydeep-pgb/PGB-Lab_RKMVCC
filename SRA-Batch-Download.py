import subprocess
import os

def read_srr_accessions(file_path):
    """
    Read SRR accessions from a file.

    Parameters:
    file_path (str): Path to the file containing SRR accession IDs.

    Returns:
    list: List of SRR accession IDs.
    """
    with open(file_path, 'r') as file:
        srr_accessions = [line.strip() for line in file if line.strip()]
    return srr_accessions

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
        print(f"Running command: {' '.join(command)}")
        try:
            subprocess.run(command, check=True)
            print(f"Successfully downloaded {srr}")
        except subprocess.CalledProcessError as e:
            print(f"Failed to download {srr}: {e}")

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
        sra_file = os.path.join(sra_dir, srr + ".sra")
        if not os.path.isfile(sra_file):
            print(f"SRA file for {srr} not found in {sra_dir}. Skipping.")
            continue
        
        command = ["fastq-dump", "--split-files", "--gzip", "--outdir", fastq_dir, sra_file]
        print(f"Running command: {' '.join(command)}")
        try:
            subprocess.run(command, check=True)
            print(f"Successfully converted {srr} to FASTQ files")
        except subprocess.CalledProcessError as e:
            print(f"Failed to convert {srr} to FASTQ files: {e}")

if __name__ == "__main__":
    # Path to the file containing SRR accession IDs
    srr_file_path = "srr_accessions.txt"
    
    # Directory to download SRA files to
    sra_dir = "/path/to/downloaded/sra/files"
    
    # Directory to store FASTQ files
    fastq_dir = "/path/to/output/fastq/files"

    # Create directories if they don't exist
    os.makedirs(sra_dir, exist_ok=True)
    os.makedirs(fastq_dir, exist_ok=True)

    # Read SRR accessions from file
    srr_accessions = read_srr_accessions(srr_file_path)

    # Download SRA files
    prefetch_sra_files(srr_accessions, sra_dir)

    # Convert SRA files to FASTQ files
    fastq_dump_paired_end(srr_accessions, sra_dir, fastq_dir)
