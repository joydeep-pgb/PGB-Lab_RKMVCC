import os
import subprocess

# Define the directory where your BAM files are located
bam_dir = "/media/joydeep/Elements/Mango_Hongyu/mapped/"

# Define the directory where you want to save the assembled transcripts
output_dir = "/media/joydeep/Elements/Mango_Hongyu/Assembled/"

# Define the path to the reference annotation GTF file
reference_gtf = "/media/joydeep/Elements/Mango_Hongyu/mango_TA4_gene_models.gff"

# Create the output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# List all the BAM files in the input directory
bam_files = [f for f in os.listdir(bam_dir) if f.endswith(".bam")]

# Loop through each BAM file and run StringTie
for bam_file in bam_files:
    sample_name = bam_file.replace(".bam", "")
    input_bam = os.path.join(bam_dir, bam_file)
    output_gtf = os.path.join(output_dir, f"{sample_name}.gtf")
    abundance_file = os.path.join(output_dir, f"{sample_name}_abundance.gtf")

    # Run StringTie command for assembly with reference annotations
    stringtie_cmd = f"stringtie -p 4 {input_bam} -o {output_gtf} -e -A {abundance_file} -G {reference_gtf} -f 0.01 -m 200 -a 10 -j 1 -c 1 -g 50 -M 1.0"
    subprocess.run(stringtie_cmd, shell=True)

    print(f"Assembled: {sample_name}")

print("Assembly complete with refference transcriptsonly.")

