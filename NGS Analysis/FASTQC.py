import os
import subprocess

# Define input and output folders
input_folder = "/mnt/WD_Blue/Sorghum/TEST/FASTQ/"
output_folder = "/mnt/WD_Blue/Sorghum/TEST/fastqc/"

# Create output folder if it doesn't exist
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Get a list of all FASTQ files in the input folder
fastq_files = [file for file in os.listdir(input_folder) if file.endswith('.fastqsanger.gz')]

# Run FastQC for each FASTQ file
for fastq_file in fastq_files:
    input_path = os.path.join(input_folder, fastq_file)
    output_path = os.path.join(output_folder, f"fastqc_{fastq_file}")

    # Run FastQC command
    command = f"fastqc -o {output_folder} {input_path}"
    subprocess.run(command, shell=True)

print("FastQC analysis completed.")
