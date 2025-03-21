import os
import subprocess

# Function to run a shell command
def run_command(command):
    try:
        subprocess.run(command, check=True, shell=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {e}")
        exit(1)

# Paths and parameters
genome_dir = 'tempstargenomedir'
genome_fasta = 'refgenome.fa'
annotation_gtf = '/data/dnb10/galaxy_db/files/c/f/4/dataset_cf44fff3-f3a9-4d8e-94ed-381a46ab71e0.dat'
num_threads = '10'

# Step 1: Generate STAR genome index
print("Generating STAR genome index...")

genome_generate_command = f"""
STAR --runMode genomeGenerate \
     --genomeDir {genome_dir} \
     --genomeFastaFiles {genome_fasta} \
     --sjdbOverhang 100 \
     --sjdbGTFfile {annotation_gtf} \
     --sjdbGTFfeatureExon exon \
     --sjdbGTFtagExonParentTranscript Parent \
     --genomeSAindexNbases 14 \
     --runThreadN {num_threads}
"""
run_command(genome_generate_command)

print("STAR genome index complete.")
