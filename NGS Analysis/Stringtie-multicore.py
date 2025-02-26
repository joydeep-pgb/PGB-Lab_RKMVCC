import os
import subprocess
import multiprocessing
import logging

def run_stringtie(bam_file, bam_dir, output_dir, reference_gtf):
    """Runs StringTie for a single BAM file."""
    sample_name = bam_file.replace(".bam", "")
    input_bam = os.path.join(bam_dir, bam_file)
    output_gtf = os.path.join(output_dir, f"{sample_name}.gtf")
    abundance_file = os.path.join(output_dir, f"{sample_name}_abundance.gtf")

    logging.info(f"Processing file: {bam_file}")

    stringtie_cmd = f"stringtie -p 1 -G {reference_gtf} -e -o {output_gtf} {input_bam} -A {abundance_file} -f 0.01 -m 200 -a 10 -j 1 -c 1 -g 50 -M 1.0"
    try:
        subprocess.run(stringtie_cmd, shell=True, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        logging.info(f"Assembled: {sample_name}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error assembling {sample_name}: {e}")
        logging.error(f"Command: {stringtie_cmd}")
        logging.error(f"Stderr: {e.stderr.decode()}")
        logging.error(f"Stdout: {e.stdout.decode()}")
    except Exception as e:
        logging.error(f"Unexpected error assembling {sample_name}: {e}")
        logging.error(f"Command: {stringtie_cmd}")

def main(num_processes=None):
    """Main function to run StringTie in parallel."""
    # Define the directory where your BAM files are located
    bam_dir = "/home/pgb-lab/Documents/Sorghum_MetaDEG/Salt/Mapped/"

    # Define the directory where you want to save the assembled transcripts
    output_dir = "/home/pgb-lab/Documents/Sorghum_MetaDEG/Salt/Assembled/"

    # Define the path to the reference annotation GTF file
    reference_gtf = "/home/pgb-lab/Documents/Sorghum_MetaDEG/GFF/Sbicolor.gff3"

    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Configure logging to both file and console
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler('stringtie_assembly.log'),
            logging.StreamHandler()
        ]
    )

    # List all the BAM files in the input directory
    bam_files = [f for f in os.listdir(bam_dir) if f.endswith(".bam")]

    if not bam_files:
        logging.warning(f"No BAM files found in {bam_dir}")
        return

    # Create a pool of worker processes
    with multiprocessing.Pool(processes=num_processes) as pool:
        # Prepare arguments for each call to run_stringtie
        tasks = [(bam_file, bam_dir, output_dir, reference_gtf) for bam_file in bam_files]
        # Run StringTie in parallel
        pool.starmap(run_stringtie, tasks)

    logging.info("Assembly process completed.")

if __name__ == "__main__":
    num_processes_to_use = 10
    main(num_processes=num_processes_to_use)
