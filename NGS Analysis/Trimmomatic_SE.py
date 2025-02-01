import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

def trimmomatic_single(input_dir, output_dir):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    
    # Create the output directory if it doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # List all the single-end FASTQ files in the input directory
    fastq_files = list(input_dir.glob("*.fastqsanger.gz"))

    with ThreadPoolExecutor(max_workers=4) as executor:
        for fastq_file in fastq_files:
            sample_name = fastq_file.stem
            output_path_paired = output_dir / f"{sample_name}_paired.fastqsanger.gz"
            output_path_unpaired = output_dir / f"{sample_name}_unpaired.fastqsanger.gz"

            # Run Trimmomatic command for single-end reads in parallel
            trim_cmd = [
                "trimmomatic", 
                "SE", "-threads", "4", str(fastq_file), 
                str(output_path_paired), str(output_path_unpaired), 
                "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:20", "MINLEN:20"
            ]
            
            executor.submit(subprocess.run, trim_cmd)
    
    print("Trimming complete.")

if __name__ == "__main__":
    input_dir = "/mnt/WD_Blue/Sorghum/TEST/FASTQ/"
    output_dir = "/mnt/WD_Blue/Sorghum/TEST/Trimmed/"

    trimmomatic_single(input_dir, output_dir)
