import os
import subprocess

def trimmomatic_paired(input_dir, output_dir):
    # Create the output and summary directories if they don't exist
    os.makedirs(output_dir, exist_ok=True)
    #os.makedirs(summary_dir, exist_ok=True)

    # List all the paired-end FASTQ files in the input directory
    fastq_files_1 = [f for f in os.listdir(input_dir) if f.endswith("_forward.fastqsanger.gz")]

    # Loop through each paired-end sample and run Trimmomatic
    for fastq_file_1 in fastq_files_1:
        sample_name = fastq_file_1.replace("_forward.fastqsanger.gz", "")
        fastq_file_2 = fastq_file_1.replace("_forward.fastqsanger.gz", "_reverse.fastqsanger.gz")

        input_path_1 = os.path.join(input_dir, fastq_file_1)
        input_path_2 = os.path.join(input_dir, fastq_file_2)
        output_path_1_paired = os.path.join(output_dir, f"{sample_name}_forward_paired.fastqsanger.gz")
        output_path_1_unpaired = os.path.join(output_dir, f"{sample_name}_forward_unpaired.fastqsanger.gz")
        output_path_2_paired = os.path.join(output_dir, f"{sample_name}_reverse_paired.fastqsanger.gz")
        output_path_2_unpaired = os.path.join(output_dir, f"{sample_name}_reverse_unpaired.fastqsanger.gz")

        # Run Trimmomatic command for paired-end reads
        trim_cmd = f"java -jar /mnt/WD_Blue/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 4 {input_path_1} {input_path_2} {output_path_1_paired} {output_path_1_unpaired} {output_path_2_paired} {output_path_2_unpaired} LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:20"
        subprocess.run(trim_cmd, shell=True)

        print(f"Processed: {sample_name}")

    print("Trimming complete.")

if __name__ == "__main__":
    input_dir = "/mnt/WD_Blue/Sorghum/TEST/FASTQ/"
    output_dir = "/mnt/WD_Blue/Sorghum/TEST/Trimmed/"

    trimmomatic_paired(input_dir, output_dir)
