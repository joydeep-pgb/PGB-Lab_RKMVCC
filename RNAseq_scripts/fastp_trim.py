import os
import subprocess

# Define the directory where your paired-end FASTQ files are located
input_dir = "/media/joydeep/Elements/Mango_Hongyu/fastq_files/"

# Define the directory where you want to save the trimmed files
output_dir = "/media/joydeep/Elements/Mango_Hongyu/trimmed/"

# Define the directory for reports (HTML and JSON)
report_dir = "/media/joydeep/Elements/Mango_Hongyu/fastp_reports/"

# Create the output and report directories if they don't exist
os.makedirs(output_dir, exist_ok=True)
os.makedirs(report_dir, exist_ok=True)

# List all the paired-end FASTQ files in the input directory
fastq_files_1 = [f for f in os.listdir(input_dir) if f.endswith("_forward.fastqsanger")]

# Loop through each paired-end sample and run fastp
for fastq_file_1 in fastq_files_1:
    sample_name = fastq_file_1.replace("_forward.fastqsanger", "")
    fastq_file_2 = fastq_file_1.replace("_forward.fastqsanger", "_reverse.fastqsanger")

    input_path_1 = os.path.join(input_dir, fastq_file_1)
    input_path_2 = os.path.join(input_dir, fastq_file_2)
    output_path_1 = os.path.join(output_dir, fastq_file_1)
    output_path_2 = os.path.join(output_dir, fastq_file_2)
    html_report = os.path.join(report_dir, f"{sample_name}_report.html")
    json_report = os.path.join(report_dir, f"{sample_name}_report.json")

    # Run fastp command for paired-end inputs with HTML and JSON reports
    cmd = f"fastp -w 4 -i {input_path_1} -I {input_path_2} -o {output_path_1} -O {output_path_2} -h {html_report} -j {json_report}"
    subprocess.run(cmd, shell=True)

    print(f"Processed: {sample_name}")

print("Processing complete.")
