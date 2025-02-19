import os
import subprocess

# Define directories
input_dir = "/mnt/d/Phaseolus_Tissue/Trimmed/"
output_dir = "/mnt/d/Phaseolus_Tissue/Mapped/"
index_path = "/mnt/d/Phaseolus_Tissue/Phaseolus_Index/Pvulgaris_Index"
summary_dir = "/mnt/d/Phaseolus_Tissue/Mapped/Summary"

# Create directories if they don't exist
os.makedirs(output_dir, exist_ok=True)
os.makedirs(summary_dir, exist_ok=True)

# List all single-end FASTQ files
fastq_files = [f for f in os.listdir(input_dir) if f.endswith(".fastq.gz")]

# Loop through samples and process
for fastq_file in fastq_files:
    sample_name = fastq_file.replace(".fastq.gz", "")
    input_path = os.path.join(input_dir, fastq_file)
    sorted_bam_output = os.path.join(output_dir, f"{sample_name}_sorted.bam")
    summary_file = os.path.join(summary_dir, f"{sample_name}_summary.txt")

    # Run HISAT2 mapping
    hisat2_cmd = f"hisat2 -p 8 --dta -x {index_path} -U {input_path} --new-summary --summary-file {summary_file}"
    print(f"Processing: {sample_name}")
    
    # Run HISAT2 and pipe output directly to samtools for sorting
    samtools_cmd = f"samtools view -@ 8 -bS - | samtools sort -@ 8 -o {sorted_bam_output} -"
    
    # Combine both commands using a pipe
    full_cmd = f"{hisat2_cmd} | {samtools_cmd}"
    subprocess.run(full_cmd, shell=True, executable="/bin/bash")  # Ensure Bash is used for piping

print("Mapping and sorting complete.")
