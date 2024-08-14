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
read1 = '/data/dnb10/galaxy_db/files/e/d/3/dataset_ed355746-3420-458d-a1c4-a449a03c5c5b.dat'
read2 = '/data/dnb10/galaxy_db/files/5/0/a/dataset_50a28709-c073-4201-af55-5430f84acc8b.dat'
output_prefix = '/path/to/output/sample'
num_threads = os.getenv('GALAXY_SLOTS', 4)

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

# Step 2: Align RNA-seq reads
print("Aligning RNA-seq reads...")
align_command = f"""
STAR --runThreadN {num_threads} \
     --genomeLoad NoSharedMemory \
     --genomeDir {genome_dir} \
     --readFilesIn {read1} {read2} \
     --readFilesCommand zcat \
     --outFileNamePrefix {output_prefix} \
     --outSAMtype BAM SortedByCoordinate \
     --twopassMode None \
     --quantMode - \
     --outSAMattrIHstart 1 \
     --outSAMattributes NH HI AS nM ch \
     --outSAMprimaryFlag OneBestScore \
     --outSAMmapqUnique 60 \
     --outSAMunmapped Within \
     --outBAMsortingThreadN {num_threads} \
     --outBAMsortingBinsN 50 \
     --winAnchorMultimapNmax 50
"""
run_command(align_command)

print("RNA-seq read alignment completed.")
