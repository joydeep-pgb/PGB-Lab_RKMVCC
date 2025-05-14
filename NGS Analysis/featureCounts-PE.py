import os
import subprocess
import glob

# Configuration
input_dir = "/media/pgb-lab/One_HDD/Pvulgaris_DEG/NEW_SRA/Drought_Root/Mapped/"      # Directory containing BAM files
output_dir = "/media/pgb-lab/One_HDD/Pvulgaris_DEG/NEW_SRA/Drought_Root/featureCounts/"        # Output directory for count files
gtf_file = "/media/pgb-lab/One_HDD/Pvulgaris_DEG/NEW_SRA/Pvulgaris_Index/Pvulgaris.gff3"  # Path to GTF annotation file
threads = 4                           # Number of threads to use
paired_end = True                     # Set to False for single-end data

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Get list of BAM files (exclude index files)
bam_paths = glob.glob(os.path.join(input_dir, "*.bam"))
bam_paths = [bam for bam in bam_paths if not bam.endswith('.bai')]

if not bam_paths:
    raise FileNotFoundError(f"No BAM files found in {input_dir}")

# Set output file path
output_file = "combined_counts.txt"
output_path = os.path.join(output_dir, output_file)

# Base command construction
cmd = [
    "featureCounts",
    "-a", gtf_file,
    "-o", output_path,
    "-T", str(threads),
    "-s", "0",                # Strandedness: 0=unstranded
    "-t", "gene",             # Feature type to count
    "-g", "ID",               # GTF attribute to use
    "-Q", "0",                # Minimum mapping quality
    "--minOverlap", "1",      # Minimum overlapping bases
    "--fracOverlap", "0",     # Minimum fractional read overlap
    "--fracOverlapFeature", "0",  # Minimum fractional feature overlap
    "-C"                      # Exclude chimeric reads
]

# Add paired-end options if needed
if paired_end:
    cmd += ["-p", "--countReadPairs"]

# Add BAM files to command
cmd += bam_paths

# Print and execute command
print("Running:", ' '.join(cmd))
try:
    result = subprocess.run(cmd, check=True, text=True, stderr=subprocess.PIPE)
    print("\nSuccess! Output saved to:", output_path)
    print("featureCounts log:")
    print(result.stderr)
except subprocess.CalledProcessError as e:
    print(f"\nError running featureCounts (code {e.returncode}):")
    print(e.stderr)
    raise
except Exception as e:
    print(f"Unexpected error: {str(e)}")
    raise

print("\nProcessing completed successfully!")