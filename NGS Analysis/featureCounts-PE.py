import os
import subprocess

# Configuration
input_dir = "/path/to/bam/files"      # Directory containing BAM files
output_dir = "/path/to/output"        # Output directory for count files
gtf_file = "/path/to/annotation.gtf"  # Path to your GTF annotation file
threads = 4                           # Number of threads to use

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Get list of BAM files in input directory
bam_files = [f for f in os.listdir(input_dir) if f.endswith('.bam')]

for bam_file in bam_files:
    # Input BAM path
    input_bam = os.path.join(input_dir, bam_file)
    
    # Output file name (using BAM file name as base)
    output_file = os.path.splitext(bam_file)[0] + "_counts.txt"
    output_path = os.path.join(output_dir, output_file)
    
    # Construct featureCounts command
    cmd = [
        "featureCounts",
        "-a", gtf_file,
        "-o", output_path,
        "-p",                          # For paired-end data
        "--countReadPairs",            # Count read pairs instead of reads
        "-B",                          # Require both ends to be mapped
        "-T", str(threads),            # Number of threads
        input_bam                      # Input BAM file
    ]
    
    # Print command (for debugging)
    print("Running:", ' '.join(cmd))
    
    # Run command
    try:
        subprocess.run(cmd, check=True)
        print(f"Successfully processed {bam_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error processing {bam_file}: {e}")

print("All files processed!")
