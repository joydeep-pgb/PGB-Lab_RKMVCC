import os
import subprocess

# Define directories
input_dir = "/mnt/d/LAB_Data/Phaseolus_DEG/Phaseolus_Leaf/SUS/Drought/Trimmed/"
output_dir = "/mnt/d/LAB_Data/Phaseolus_DEG/Phaseolus_Leaf/SUS/Drought/Mapped/"
index_path = "/mnt/d/LAB_Data/Phaseolus_Tissue/Phaseolus_Index/Pvulgaris_Index"
summary_dir = "/mnt/d/LAB_Data/Phaseolus_DEG/Phaseolus_Leaf/SUS/Drought/Mapped/Summary/"

# Create directories if they don't exist
os.makedirs(output_dir, exist_ok=True)
os.makedirs(summary_dir, exist_ok=True)

# List all forward reads
fastq_files_1 = [f for f in os.listdir(input_dir) if f.endswith("_1.fastq.gz")]

# Loop through paired-end samples and process them
for fastq_file_1 in fastq_files_1:
    sample_name = fastq_file_1.replace("_1.fastq.gz", "")
    fastq_file_2 = fastq_file_1.replace("_1.fastq.gz", "_2.fastq.gz")

    input_path_1 = os.path.join(input_dir, fastq_file_1)
    input_path_2 = os.path.join(input_dir, fastq_file_2)
    sorted_bam_output = os.path.join(output_dir, f"{sample_name}_sorted.bam")
    summary_file = os.path.join(summary_dir, f"{sample_name}_summary.txt")

    # Run HISAT2, pipe output to samtools for BAM conversion & sorting
    hisat2_cmd = f"hisat2 -p 8 --dta -x {index_path} -1 {input_path_1} -2 {input_path_2} --new-summary --summary-file {summary_file}"
    samtools_cmd = f"samtools view -@ 8 -bS - | samtools sort -@ 8 -o {sorted_bam_output} -"

    full_cmd = f"{hisat2_cmd} | {samtools_cmd}"
    
    print(f"Processing: {sample_name}")
    subprocess.run(full_cmd, shell=True, executable="/bin/bash")

    # Index the sorted BAM file
    samtools_index_cmd = f"samtools index {sorted_bam_output}"
    subprocess.run(samtools_index_cmd, shell=True)

print("Mapping, sorting, and indexing complete.")
