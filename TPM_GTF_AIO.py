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
  with open(file_path, 'r') as file:
    lines = file.readlines()

  data = [line.strip().split('\t') for line in lines if 'TPM' in line]

  result = []
  for line in data:
    transcript_id = re.search(r'transcript_id "([^"]*)"', line[8]).group(1)
    TPM_value = re.search(r'TPM "([^"]*)"', line[8]).group(1)
    result.append((transcript_id, TPM_value))

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
      with open(output_file, 'w') as out_file:
        out_file.write(f"transcript_id\tTPM\n")
        for row in data:
          out_file.write(f"{row[0]}\t{row[1]}\n")

if __name__ == "__main__":
  if len(sys.argv) < 3:
    print("Usage: python script.py gtf_folder output_folder")
    sys.exit(1)

  gtf_folder = sys.argv[1]
  output_folder = sys.argv[2]

  extract_from_folder(gtf_folder, output_folder)
