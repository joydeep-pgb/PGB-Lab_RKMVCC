import sys
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
                        transcript_id = re.search(r'transcript_id "([^"]*)"', attributes).group(1)
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
                    out_file.write(f"transcript_id\t{filename[:-4]}\n")  # Header changed to the file name
                    for row in data:
                        out_file.write(f"{row[0]}\t{row[1]}\n")
            except Exception as e:
                print(f"An error occurred while writing to file {output_file}: {e}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script.py gtf_folder output_folder")
        sys.exit(1)

    gtf_folder = sys.argv[1]
    output_folder = sys.argv[2]

    extract_from_folder(gtf_folder, output_folder)
