import os
import subprocess
import multiprocessing
from datetime import datetime

def trimmomatic_paired_sample(input_path_1, input_path_2, output_path_1_paired, output_path_1_unpaired, output_path_2_paired, output_path_2_unpaired, log_path):
    """Trims a single paired-end sample with logging."""
    sample_name = os.path.basename(input_path_1).replace("_1.fastq.gz", "")

    def log(message):
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        log_message = f"[{timestamp}] {message}"
        print(log_message)
        with open(log_path, "a") as f:
            f.write(log_message + "\n")

    try:
        log(f"Starting processing for {sample_name}")

        cmd_list = [
            "trimmomatic", "PE", "-threads", "1",
            input_path_1, input_path_2,
            output_path_1_paired, output_path_1_unpaired,
            output_path_2_paired, output_path_2_unpaired,
            "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:20", "MINLEN:20"
        ]

        log(f"Executing: {' '.join(cmd_list)}")

        with open(log_path, "a") as f:
            f.write("\n=== TRIMMOMATIC OUTPUT ===\n")
            subprocess.run(
                cmd_list,
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

def trimmomatic_paired_parallel(input_dir, output_dir, num_processes):
    """Trims multiple paired-end samples in parallel with logging."""
    input_dir = os.path.abspath(input_dir)
    output_dir = os.path.abspath(output_dir)
    log_dir = os.path.join(output_dir, "process_logs")
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)

    fastq_files_1 = [f for f in os.listdir(input_dir) if f.endswith("_1.fastq.gz")]
    samples = []

    for fastq_file_1 in fastq_files_1:
        sample_name = fastq_file_1.replace("_1.fastq.gz", "")
        fastq_file_2 = fastq_file_1.replace("_1.fastq.gz", "_2.fastq.gz")

        input_path_1 = os.path.join(input_dir, fastq_file_1)
        input_path_2 = os.path.join(input_dir, fastq_file_2)
        output_path_1_paired = os.path.join(output_dir, f"{sample_name}_1_paired.fastq.gz")
        output_path_1_unpaired = os.path.join(output_dir, f"{sample_name}_1_unpaired.fastq.gz")
        output_path_2_paired = os.path.join(output_dir, f"{sample_name}_2_paired.fastq.gz")
        output_path_2_unpaired = os.path.join(output_dir, f"{sample_name}_2_unpaired.fastq.gz")
        log_path = os.path.join(log_dir, f"{sample_name}.log")

        samples.append((
            input_path_1, input_path_2,
            output_path_1_paired, output_path_1_unpaired,
            output_path_2_paired, output_path_2_unpaired,
            log_path
        ))

    with multiprocessing.Pool(processes=num_processes) as pool:
        results = []
        for sample in samples:
            results.append(pool.apply_async(trimmomatic_paired_sample, sample))

        success = True
        for result in results:
            if not result.get():
                success = False

        if success:
            print("All samples processed successfully.")
        else:
            print("Some samples encountered errors during processing. Check the log files.")

if __name__ == "__main__":
    input_dir = "/media/pgb-lab/One_HDD/Pvulgaris_DEG/NEW_SRA/Drought_Root/FASTQ/"
    output_dir = "/media/pgb-lab/One_HDD/Pvulgaris_DEG/NEW_SRA/Drought_Root/Trimmed/NEW/"
    num_processes = 4

    trimmomatic_paired_parallel(input_dir, output_dir, num_processes)