import os
import subprocess
import multiprocessing

def trimmomatic_single_sample(input_path, output_path):
    """Trims a single FASTQ file."""
    sample_name = os.path.basename(input_path).replace(".fastq.gz", "").replace(".fq.gz", "")
    trim_cmd = [
        "trimmomatic", "SE",
        "-threads", "1",
        input_path,
        output_path,
        "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:20", "MINLEN:20"
    ]
    try:
        subprocess.run(trim_cmd, check=True)
        print(f"Processed: {sample_name}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error processing {sample_name}: {e}")
        return False
    except Exception as e:
        print(f"Unexpected error processing {sample_name}: {e}")
        return False

def trimmomatic_parallel(input_dir, output_dir, num_processes):
    """Trims multiple FASTQ files in parallel."""

    input_dir = os.path.abspath(input_dir)  # Use absolute paths (recommended)
    output_dir = os.path.abspath(output_dir) # Use absolute paths (recommended)
    os.makedirs(output_dir, exist_ok=True)

    fastq_files = [f for f in os.listdir(input_dir) if f.endswith((".fastq.gz", ".fq.gz"))]
    file_paths = [os.path.join(input_dir, fastq_file) for fastq_file in fastq_files]

    with multiprocessing.Pool(processes=num_processes) as pool:
        results = []
        for file_path in file_paths:
            sample_name = os.path.basename(file_path).replace(".fastq.gz", "").replace(".fq.gz", "")
            output_path = os.path.join(output_dir, f"{sample_name}_trimmed.fastq.gz")
            results.append(pool.apply_async(trimmomatic_single_sample, (file_path, output_path)))

        success = True
        for result in results:
            if not result.get():
                success = False

        if success:
            print("Trimming complete.")
        else:
            print("Trimming completed with some errors. Check the logs.")

if __name__ == "__main__":
    input_dir = "/mnt/Blue_Drive/NGS_WD/Sorghum_MetaDEG/FASTQ/"
    output_dir = "/mnt/Blue_Drive/NGS_WD/Sorghum_MetaDEG/FASTQ_trimmed/"
    num_processes = multiprocessing.cpu_count()  # Use all available cores
    num_processes = 4

    trimmomatic_parallel(input_dir, output_dir, num_processes)
