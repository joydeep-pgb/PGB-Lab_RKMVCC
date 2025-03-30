import os
import subprocess
import shutil
import multiprocessing
import shlex
from datetime import datetime

# Configuration
input_dir = "/home/pgb-lab/Documents/Sorghum_MetaDEG/Salt/Trimmed/"
output_dir = "/media/pgb-lab/One Touch/Sorghum_MetaDEG/Salt/STAR_Mapped/"
index_path = "/home/pgb-lab/Documents/Sorghum_MetaDEG/STAR_Genome/"
summary_dir = os.path.join(output_dir, "Summary")
log_dir = os.path.join(summary_dir, "process_logs")

# Resource allocation
star_threads = 4
max_processes = 4  # Uses 16 cores (8 threads × 2 processes)

# Create directories
os.makedirs(output_dir, exist_ok=True)
os.makedirs(summary_dir, exist_ok=True)
os.makedirs(log_dir, exist_ok=True)

# Get sample list (SINGLE-END MODIFICATION HERE)
fastq_files = [f for f in os.listdir(input_dir) if f.endswith(".fastq.gz")]
samples = []
for fq in fastq_files:
    sample = fq.replace(".fastq.gz", "")
    samples.append((
        sample,
        os.path.join(input_dir, fq)
    ))

def process_sample(sample):
    """Process a single sample with logging and error handling"""
    sample_name, fq_path = sample  # SINGLE-END MODIFICATION HERE
    log_path = os.path.join(log_dir, f"{sample_name}.log")
    sorted_bam = os.path.join(output_dir, f"{sample_name}_sorted.bam")
    summary_file = os.path.join(summary_dir, f"{sample_name}_summary.txt")

    def log(message):
        """Helper function for logging"""
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        log_message = f"[{timestamp}] {message}"
        print(log_message)
        with open(log_path, "a") as f:
            f.write(log_message + "\n")

    try:
        log(f"Starting processing for {sample_name}")
        
        # STAR Command (SINGLE-END MODIFICATION HERE)
        star_cmd = (
            f"STAR --runThreadN {star_threads} "
            f"--genomeDir {shlex.quote(index_path)} "
            f"--readFilesIn {shlex.quote(fq_path)} "  # Single file input
            f"--readFilesCommand zcat "
            f"--outSAMtype BAM SortedByCoordinate "
            f"--outTmpDir /tmp/{sample_name}_STARtmp "
            f"--outFileNamePrefix {shlex.quote(os.path.join(output_dir, sample_name + '_'))}"
        )

        log(f"Executing: {star_cmd}")
        with open(log_path, "a") as f:
            f.write("\n=== STAR OUTPUT ===\n")
            subprocess.run(
                star_cmd,
                shell=True,
                check=True,
                stdout=f,
                stderr=subprocess.STDOUT,
                text=True
            )

        # Handle STAR outputs
        star_bam = os.path.join(output_dir, f"{sample_name}_Aligned.sortedByCoord.out.bam")
        if os.path.exists(star_bam):
            os.rename(star_bam, sorted_bam)
            log(f"Renamed BAM file to {sorted_bam}")
        else:
            raise FileNotFoundError(f"STAR output BAM missing: {star_bam}")

        # Move summary file
        star_log = os.path.join(output_dir, f"{sample_name}_Log.final.out")
        if os.path.exists(star_log):
            shutil.move(star_log, summary_file)
            log(f"Moved summary file to {summary_file}")
        else:
            log("Warning: Final STAR log file not found")

        # Index BAM
        log("Indexing BAM file")
        subprocess.run(
            f"samtools index {shlex.quote(sorted_bam)}",
            shell=True,
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE
        )
        log("Indexing completed successfully")

    except subprocess.CalledProcessError as e:
        log(f"Command failed with error {e.returncode}: {e.stderr}")
    except Exception as e:
        log(f"Critical error: {str(e)}")
    finally:
        log(f"Processing completed for {sample_name}")

if __name__ == "__main__":
    print(f"Starting parallel processing of {len(samples)} samples")
    with multiprocessing.Pool(processes=max_processes) as pool:
        pool.map(process_sample, samples)
    print("All samples processed")