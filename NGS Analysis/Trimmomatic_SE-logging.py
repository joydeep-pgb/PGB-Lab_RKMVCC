import os
import subprocess
import multiprocessing
from datetime import datetime

def trimmomatic_single_sample(input_path, output_path, log_path):
    """Trims a single FASTQ file with logging."""
    sample_name = os.path.basename(input_path).replace(".fastq.gz", "").replace(".fq.gz", "")

    def log(message):
        """Helper function for logging"""
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        log_message = f"[{timestamp}] {message}"
        print(log_message)  # Print to terminal
        with open(log_path, "a") as f:
            f.write(log_message + "\n")

    try:
        log(f"Starting processing for {sample_name}")
        
        trim_cmd = [
            "trimmomatic", "SE",
            "-threads", "1",
            input_path,
            output_path,
            "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:20", "MINLEN:20"
        ]

        log(f"Executing: {' '.join(trim_cmd)}")

        with open(log_path, "a") as f:
            f.write("\n=== TRIMMOMATIC OUTPUT ===\n")
            subprocess.run(
                trim_cmd,
                stdout=f,
                stderr=subprocess.STDOUT,
                check=True,
                text=True
            )

        log(f"Successfully processed {sample_name}")
        return True

    except subprocess.CalledProcessError as e:
        log(f"Command failed with error {e.returncode}: {e}")
        return False
    except Exception as e:
        log(f"Critical error: {str(e)}")
        return False

def trimmomatic_parallel(input_dir, output_dir, num_processes):
    """Trims multiple FASTQ files in parallel with logging."""
    input_dir = os.path.abspath(input_dir)
    output_dir = os.path.abspath(output_dir)
    log_dir = os.path.join(output_dir, "process_logs")
    
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)

    fastq_files = [f for f in os.listdir(input_dir) if f.endswith((".fastq.gz", ".fq.gz"))]
    
    samples = []
    for fastq_file in fastq_files:
        input_path = os.path.join(input_dir, fastq_file)
        sample_name = os.path.basename(input_path).replace(".fastq.gz", "").replace(".fq.gz", "")
        output_path = os.path.join(output_dir, f"{sample_name}_trimmed.fastq.gz")
        log_path = os.path.join(log_dir, f"{sample_name}.log")
        samples.append((input_path, output_path, log_path))

    with multiprocessing.Pool(processes=num_processes) as pool:
        results = []
        for sample in samples:
            results.append(pool.apply_async(trimmomatic_single_sample, sample))

        success = True
        for result in results:
            if not result.get():
                success = False

        if success:
            print("All samples processed successfully.")
        else:
            print("Some samples encountered errors during processing. Check the log files.")

if __name__ == "__main__":
    input_dir = "/mnt/Blue_Drive/NGS_WD/Sorghum_MetaDEG/FASTQ/"
    output_dir = "/mnt/Blue_Drive/NGS_WD/Sorghum_MetaDEG/FASTQ_trimmed/"
    num_processes = 4

    trimmomatic_parallel(input_dir, output_dir, num_processes)