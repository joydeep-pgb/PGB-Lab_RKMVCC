import os
import subprocess

def trimmomatic_paired(input_dir, output_dir):
    # Create the output and summary directories if they don't exist
    os.makedirs(output_dir, exist_ok=True)
    #os.makedirs(summary_dir, exist_ok=True)

    # List all the paired-end FASTQ files in the input directory
    fastq_files_1 = [f for f in os.listdir(input_dir) if f.endswith("_1.fastq.gz")]

    # Loop through each paired-end sample and run Trimmomatic
    for fastq_file_1 in fastq_files_1:
        sample_name = fastq_file_1.replace("_1.fastq.gz", "")
        fastq_file_2 = fastq_file_1.replace("_1.fastq.gz", "_2.fastq.gz")

        input_path_1 = os.path.join(input_dir, fastq_file_1)
        input_path_2 = os.path.join(input_dir, fastq_file_2)
        output_path_1_paired = os.path.join(output_dir, f"{sample_name}_1_paired.fastq.gz")
        output_path_1_unpaired = os.path.join(output_dir, f"{sample_name}_1_unpaired.fastq.gz")
        output_path_2_paired = os.path.join(output_dir, f"{sample_name}_2_paired.fastq.gz")
        output_path_2_unpaired = os.path.join(output_dir, f"{sample_name}_2_unpaired.fastq.gz")

        # Run Trimmomatic command for paired-end reads
        trim_cmd = f"trimmomatic PE -threads 8 {input_path_1} {input_path_2} {output_path_1_paired} {output_path_1_unpaired} {output_path_2_paired} {output_path_2_unpaired} LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:20"
        subprocess.run(trim_cmd, shell=True)

        print(f"Processed: {sample_name}")

    print("Trimming complete.")

if __name__ == "__main__":
    input_dir = "/mnt/d/LAB_Data/Phaseolus_DEG/Phaseolus_Leaf/SUS/Cold/FASTQ/"
    output_dir = "/mnt/d/LAB_Data/Phaseolus_DEG/Phaseolus_Leaf/SUS/Cold/Trimmed/"

    trimmomatic_paired(input_dir, output_dir)
