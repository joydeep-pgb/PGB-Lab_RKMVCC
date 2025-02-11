import os
import subprocess

def trimmomatic_single(input_dir, output_dir):
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # List all single-end FASTQ files in the input directory
    fastq_files = [f for f in os.listdir(input_dir) if f.endswith(".fastq.gz") or f.endswith(".fq.gz")]

    # Loop through each single-end sample and run Trimmomatic
    for fastq_file in fastq_files:
        sample_name = fastq_file.replace(".fastq.gz", "").replace(".fq.gz", "")
        input_path = os.path.join(input_dir, fastq_file)
        output_path = os.path.join(output_dir, f"{sample_name}_trimmed.fastq.gz")

        # Run Trimmomatic command for single-end reads
        trim_cmd = f"trimmomatic SE -threads 4 {input_path} {output_path} LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:20"
        subprocess.run(trim_cmd, shell=True)

        print(f"Processed: {sample_name}")

    print("Trimming complete.")

if __name__ == "__main__":
    input_dir = "/mnt/WD_Blue/Sorghum/TEST/FASTQ/"
    output_dir = "/mnt/WD_Blue/Sorghum/TEST/Trimmed/"

    trimmomatic_single(input_dir, output_dir)
