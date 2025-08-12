import re
import os

def process_gtf(file_path):
    """
    Extracts transcript IDs and TPM values from a GTF file.

    Args:
        file_path: Path to the GTF file.

    Returns:
        A list of tuples containing (transcript_id, TPM_value).
    """
    result = []
    try:
        with open(file_path, 'r') as file:
            for line in file:
                if 'TPM' in line:
                    fields = line.strip().split('\t')
                    attributes = fields[8]
                    try:
                        transcript_id = re.search(r'gene_id "([^"]*)"', attributes).group(1)
                        TPM_value = re.search(r'TPM "([^"]*)"', attributes).group(1)
                        result.append((transcript_id, TPM_value))
                    except AttributeError:
                        print(f"Error parsing line in file {file_path}: {line.strip()}")
    except FileNotFoundError:
        print(f"File not found: {file_path}")
    except Exception as e:
        print(f"An error occurred while processing file {file_path}: {e}")

    return result

def extract_from_folder(folder_path, output_folder):
    """
    Extracts TPM values from all GTF files in a folder.

    Args:
        folder_path: Path to the folder containing GTF files.
        output_folder: Path to the folder for storing output files.
    """
    os.makedirs(output_folder, exist_ok=True)  # Create output folder if it doesn't exist
    for filename in os.listdir(folder_path):
        if filename.endswith(".gtf"):
            gtf_file = os.path.join(folder_path, filename)
            output_file = os.path.join(output_folder, f"{filename[:-4]}.txt")  # Remove .gtf extension and add .txt
            data = process_gtf(gtf_file)
            try:
                with open(output_file, 'w') as out_file:
                    out_file.write(f"gene_id\t{filename[:-4]}\n")  # Header changed to the file name
                    for row in data:
                        out_file.write(f"{row[0]}\t{row[1]}\n")
            except Exception as e:
                print(f"An error occurred while writing to file {output_file}: {e}")

if __name__ == "__main__":
    # Set your input and output folder paths here
    gtf_folder = "F:\\ABC Transporters\\Phaseolus_NEW_DEG\\NEW_SRA\\Salt_Root\\Assembled"  # Change this to your GTF files directory
    output_folder = "F:\\ABC Transporters\\Phaseolus_NEW_DEG\\NEW_SRA\\Salt_Root\\Assembled\\TPM"  # Change this to your desired output directory

    # Check if the input folder exists
    if not os.path.exists(gtf_folder):
        print(f"Input folder does not exist: {gtf_folder}")
        exit(1)
    
    print(f"Processing GTF files from: {gtf_folder}")
    print(f"Output will be saved to: {output_folder}")
    
    extract_from_folder(gtf_folder, output_folder)
    print("Processing complete!")