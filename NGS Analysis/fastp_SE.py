import os
import subprocess

# Define the directory where your FASTQ files are located
input_dir = "/media/pgblab/Elements/RNAseq_Analysis/Sorghum_lncRNA/fastq_files/"

# Define the directory where you want to save the trimmed files
output_dir = "/media/pgblab/Elements/RNAseq_Analysis/Sorghum_lncRNA/Trimmed/"

# Define the directory for reports (HTML and JSON)
report_dir = "/media/pgblab/Elements/RNAseq_Analysis/Sorghum_lncRNA/Trimmed/fastp_reports/"

# Create the output and report directories if they don't exist
os.makedirs(output_dir, exist_ok=True)
os.makedirs(report_dir, exist_ok=True)

# List all the FASTQ files in the input directory
fastq_files = [f for f in os.listdir(input_dir) if f.endswith(".fastqsanger.gz")]

# Loop through each FASTQ file and run fastp to trim it
for fastq_file in fastq_files:
    input_path = os.path.join(input_dir, fastq_file)
    output_path = os.path.join(output_dir, fastq_file.replace(".fastqsanger.gz", "_trimmed.fastqsanger.gz"))
    html_report = os.path.join(report_dir, fastq_file.replace(".fastqsanger.gz", "_report.html"))
    json_report = os.path.join(report_dir, fastq_file.replace(".fastqsanger.gz", "_report.json"))

    # Run fastp command with -h (HTML report) and -j (JSON report) options
    cmd = f"fastp -w 4 -i {input_path} -o {output_path} -h {html_report} -j {json_report}"
    subprocess.run(cmd, shell=True)

    print(f"Processed: {fastq_file}")

print("Processing complete.")
