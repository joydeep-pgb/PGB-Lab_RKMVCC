import os
import subprocess
import logging
from multiprocessing import Pool

# Configure logging
logging.basicConfig(
    filename='fastqc_analysis.log', 
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# Define input and output folders
input_folder = "/mnt/WD_Blue/Sorghum/TEST/FASTQ/"
output_folder = "/mnt/WD_Blue/Sorghum/TEST/fastqc/"

# Create output folder if it doesn't exist
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Get a list of all FASTQ files in the input folder
fastq_files = [file for file in os.listdir(input_folder) if file.endswith('.fastqsanger.gz')]

def run_fastqc(fastq_file):
    input_path = os.path.join(input_folder, fastq_file)
    command = f"fastqc -o {output_folder} {input_path}"
    
    logging.info(f"Processing file: {fastq_file}")
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    
    if result.returncode == 0:
        logging.info(f"Successfully processed: {fastq_file}")
    else:
        logging.error(f"Error processing {fastq_file}: {result.stderr}")

if __name__ == "__main__":
    num_workers = 4  # Set the number of parallel processes
    
    with Pool(num_workers) as pool:
        pool.map(run_fastqc, fastq_files)
    
    logging.info("FastQC analysis completed.")
    print("FastQC analysis completed.")
