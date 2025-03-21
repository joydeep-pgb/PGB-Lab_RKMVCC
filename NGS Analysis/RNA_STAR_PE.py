import os
import subprocess
import shutil

# Define directories
input_dir = "/home/pgb-lab/Documents/Sorghum_MetaDEG/STAR_FASTQ/"
output_dir = "/home/pgb-lab/Documents/Sorghum_MetaDEG/STAR_FASTQ/Mapped/"
index_path = "/home/pgb-lab/Documents/Sorghum_MetaDEG/STAR_Genome/Index/"  # STAR genome index directory
summary_dir = "/home/pgb-lab/Documents/Sorghum_MetaDEG/STAR_FASTQ/Mapped/Summary/"

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

    # Construct STAR command
    star_cmd = (
        f"STAR --runThreadN 10 "
        f"--genomeDir {index_path} "
        f"--readFilesIn {input_path_1} {input_path_2} "
        f"--readFilesCommand zcat "
        f"--outSAMtype BAM SortedByCoordinate "
        f"--outFileNamePrefix {os.path.join(output_dir, f'{sample_name}_')}"
    )

    print(f"Processing: {sample_name}")
    try:
        subprocess.run(star_cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running STAR for {sample_name}: {e}")
        continue

    # Handle output files
    try:
        # Rename sorted BAM file
        star_bam = os.path.join(output_dir, f"{sample_name}_Aligned.sortedByCoord.out.bam")
        if os.path.exists(star_bam):
            os.rename(star_bam, sorted_bam_output)
        else:
            print(f"Error: STAR output BAM not found for {sample_name}")
            continue

        # Move summary file
        star_log = os.path.join(output_dir, f"{sample_name}_Log.final.out")
        if os.path.exists(star_log):
            shutil.move(star_log, summary_file)
        else:
            print(f"Warning: STAR log file not found for {sample_name}")

        # Index sorted BAM
        subprocess.run(f"samtools index {sorted_bam_output}", shell=True, check=True)
    except Exception as e:
        print(f"Error processing output files for {sample_name}: {e}")
        continue

print("STAR mapping, sorting, and indexing complete.")
