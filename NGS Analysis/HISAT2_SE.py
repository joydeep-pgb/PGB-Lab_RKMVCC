import os
import subprocess

# Define the directory where your single-end FASTQ files are located
input_dir = "/media/joydeep/Elements/RNAseq_Analysis/Sorghum_lncRNA/Trimmomatic/"

# Define the directory where you want to save the mapped BAM files
output_dir = "/media/joydeep/Elements/RNAseq_Analysis/Sorghum_lncRNA/Mapped/"

# Define the location of the HISAT2 index
index_path = "/media/joydeep/Elements/RNAseq_Analysis/Sorghum_lncRNA/HISAT2_Index/sbi_genome"

# Define the directory for HISAT2 summary files
summary_dir = "/media/joydeep/Elements/RNAseq_Analysis/Sorghum_lncRNA/Mapped/Summary/"

# Create the output and summary directories if they don't exist
os.makedirs(output_dir, exist_ok=True)
os.makedirs(summary_dir, exist_ok=True)

# List all the single-end FASTQ files in the input directory
fastq_files = [f for f in os.listdir(input_dir) if f.endswith(".fastqsanger.gz")]

# Loop through each single-end sample and run HISAT2
for fastq_file in fastq_files:
    sample_name = fastq_file.replace(".fastqsanger.gz", "")
    input_path = os.path.join(input_dir, fastq_file)
    output_sam = os.path.join(output_dir, f"{sample_name}.sam")
    summary_file = os.path.join(summary_dir, f"{sample_name}_summary.txt")

    # Run HISAT2 command for single-end inputs with --new-summary and --summary-file options
    hisat2_cmd = f"hisat2 -p 6 --dta -x {index_path} -U {input_path} -S {output_sam} --new-summary --summary-file {summary_file}"
    subprocess.run(hisat2_cmd, shell=True)

    print(f"Processed: {sample_name}")

    # Convert SAM to BAM
    bam_output = os.path.join(output_dir, f"{sample_name}.bam")
    samtools_view_cmd = f"samtools view -bS {output_sam} > {bam_output}"
    subprocess.run(samtools_view_cmd, shell=True)

    # Remove the original SAM file
    os.remove(output_sam)

    # Sort BAM file
    sorted_bam_output = os.path.join(output_dir, f"{sample_name}_sorted.bam")
    samtools_sort_cmd = f"samtools sort -o {sorted_bam_output} {bam_output}"
    subprocess.run(samtools_sort_cmd, shell=True)

    # Remove the unsorted BAM file
    os.remove(bam_output)

print("Mapping and sorting complete.")
