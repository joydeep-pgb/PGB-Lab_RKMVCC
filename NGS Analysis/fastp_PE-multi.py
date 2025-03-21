import os
import subprocess
import logging
from multiprocessing import Pool

# Configure logging
logging.basicConfig(
    filename='fastp_analysis.log', 
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# Define directories
input_dir = "/run/media/joydeep/Elements/Laboratory_Data_Analysis/Mango_Metaanalysis_lncRNA/SRA-Files/Green_fruit/"
output_dir = "/run/media/joydeep/Elements/Laboratory_Data_Analysis/Mango_Metaanalysis_lncRNA/fastp_reslut/SRA/"
report_dir = "/run/media/joydeep/Elements/Laboratory_Data_Analysis/Mango_Metaanalysis_lncRNA/fastp_reslut/fastp-reports/"

# Create necessary directories
os.makedirs(output_dir, exist_ok=True)
os.makedirs(report_dir, exist_ok=True)

# Get list of paired-end FASTQ files
fastq_files_1 = [f for f in os.listdir(input_dir) if f.endswith("_forward.fastqsanger.gz")]

def run_fastp(fastq_file_1):
    sample_name = fastq_file_1.replace("_forward.fastqsanger.gz", "")
    fastq_file_2 = fastq_file_1.replace("_forward.fastqsanger.gz", "_reverse.fastqsanger.gz")
    
    input_path_1 = os.path.join(input_dir, fastq_file_1)
    input_path_2 = os.path.join(input_dir, fastq_file_2)
    output_path_1 = os.path.join(output_dir, f"{sample_name}_trimmed_forward.fastqsanger.gz")
    output_path_2 = os.path.join(output_dir, f"{sample_name}_trimmed_reverse.fastqsanger.gz")
    html_report = os.path.join(report_dir, f"{sample_name}_report.html")
    json_report = os.path.join(report_dir, f"{sample_name}_report.json")
    
    # Construct fastp command
    fastp_cmd = (
        f"fastp -i {input_path_1} -I {input_path_2} -o {output_path_1} -O {output_path_2} "
        f"-h {html_report} -j {json_report} --thread 16 -R 'Fastp report for {sample_name}'"
    )
    
    logging.info(f"Processing: {sample_name}")
    result = subprocess.run(fastp_cmd, shell=True, capture_output=True, text=True)
    
    if result.returncode == 0:
        logging.info(f"Successfully processed: {sample_name}")
    else:
        logging.error(f"Error processing {sample_name}: {result.stderr}")

if __name__ == "__main__":
    num_workers = 4  # Set the number of parallel processes
    
    with Pool(num_workers) as pool:
        pool.map(run_fastp, fastq_files_1)
    
    logging.info("Trimming complete.")
    print("Trimming complete.")
