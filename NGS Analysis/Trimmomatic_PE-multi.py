import os
import subprocess
import multiprocessing

def trimmomatic_paired_sample(input_path_1, input_path_2, output_path_1_paired, output_path_1_unpaired, output_path_2_paired, output_path_2_unpaired):
    """Trims a single paired-end sample."""

    trim_cmd = f"trimmomatic PE -threads 1 {input_path_1} {input_path_2} {output_path_1_paired} {output_path_1_unpaired} {output_path_2_paired} {output_path_2_unpaired} LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:20"

    try:
        subprocess.run(trim_cmd, shell=True, check=True)
        sample_name = os.path.basename(input_path_1).replace("_1.fastq.gz", "")
        print(f"Processed: {sample_name}")
        return True
    except subprocess.CalledProcessError as e:
        sample_name = os.path.basename(input_path_1).replace("_1.fastq.gz", "")
        print(f"Error processing {sample_name}: {e}")
        return False
    except Exception as e:
        sample_name = os.path.basename(input_path_1).replace("_1.fastq.gz", "")
        print(f"Unexpected error processing {sample_name}: {e}")
        return False


def trimmomatic_paired_parallel(input_dir, output_dir, num_processes):
    """Trims multiple paired-end samples in parallel."""
    input_dir = os.path.abspath(input_dir)
    output_dir = os.path.abspath(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    fastq_files_1 = [f for f in os.listdir(input_dir) if f.endswith("_1.fastq.gz")]

    with multiprocessing.Pool(processes=num_processes) as pool:
        results = []
        for fastq_file_1 in fastq_files_1:
            sample_name = fastq_file_1.replace("_1.fastq.gz", "")
            fastq_file_2 = fastq_file_1.replace("_1.fastq.gz", "_2.fastq.gz")

            input_path_1 = os.path.join(input_dir, fastq_file_1)
            input_path_2 = os.path.join(input_dir, fastq_file_2)
            output_path_1_paired = os.path.join(output_dir, f"{sample_name}_1_paired.fastq.gz")
            output_path_1_unpaired = os.path.join(output_dir, f"{sample_name}_1_unpaired.fastq.gz")
            output_path_2_paired = os.path.join(output_dir, f"{sample_name}_2_paired.fastq.gz")
            output_path_2_unpaired = os.path.join(output_dir, f"{sample_name}_2_unpaired.fastq.gz")

            results.append(pool.apply_async(trimmomatic_paired_sample, (input_path_1, input_path_2, output_path_1_paired, output_path_1_unpaired, output_path_2_paired, output_path_2_unpaired)))

        success = True
        for result in results:
            if not result.get():
                success = False

        if success:
            print("Trimming complete.")
        else:
            print("Trimming completed with some errors. Check the logs.")


if __name__ == "__main__":
    input_dir = "/mnt/Blue_Drive/NGS_WD/Sorghum_MetaDEG/Cold/Test/"  # Replace with your input directory
    output_dir = "/mnt/Blue_Drive/NGS_WD/Sorghum_MetaDEG/Cold/Test/Trimmed/" # Replace with your output directory
    num_processes = multiprocessing.cpu_count()  # Or a specific number like 4 for testing
    # num_processes = 4 # For testing

    trimmomatic_paired_parallel(input_dir, output_dir, num_processes)