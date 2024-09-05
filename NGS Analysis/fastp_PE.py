import os
import subprocess

# Define the directory where your paired-end FASTQ files are located
input_dir = "/run/media/joydeep/Elements/Laboratory_Data_Analysis/Mango_Metaanalysis_lncRNA/SRA-Files/Green_fruit/"

# Define the directory to save trimmed FASTQ files
output_dir = "/run/media/joydeep/Elements/Laboratory_Data_Analysis/Mango_Metaanalysis_lncRNA/fastp_reslut/SRA/"

# Define the directory to save fastp reports
report_dir = "/run/media/joydeep/Elements/Laboratory_Data_Analysis/Mango_Metaanalysis_lncRNA/fastp_reslut/fastp-reports/"

# Create the output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)
os.makedirs(report_dir, exist_ok=True)

# List all the paired-end FASTQ files in the input directory
fastq_files_1 = [f for f in os.listdir(input_dir) if f.endswith("_forward.fastqsanger.gz")]

# Loop through each pair of FASTQ files and run fastp
for fastq_file_1 in fastq_files_1:
    sample_name = fastq_file_1.replace("_forward.fastqsanger.gz", "")
    fastq_file_2 = fastq_file_1.replace("_forward.fastqsanger.gz", "_reverse.fastqsanger.gz")

    input_path_1 = os.path.join(input_dir, fastq_file_1)
    input_path_2 = os.path.join(input_dir, fastq_file_2)
    output_path_1 = os.path.join(output_dir, f"{sample_name}_trimmed_forward.fastqsanger.gz")
    output_path_2 = os.path.join(output_dir, f"{sample_name}_trimmed_reverse.fastqsanger.gz")
    html_report = os.path.join(report_dir, f"{sample_name}_report.html")
    json_report = os.path.join(report_dir, f"{sample_name}_report.json")

    # Run fastp command with -h, -j, and -R options for paired-end data
    fastp_cmd = f"fastp -i {input_path_1} -I {input_path_2} -o {output_path_1} -O {output_path_2} -h {html_report} -j {json_report} --thread {16} -R 'Fastp report for {sample_name}'"
    subprocess.run(fastp_cmd, shell=True)

    print(f"Processed: {sample_name}")

print("Trimming complete.")

