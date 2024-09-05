import sys
import re
import os

def process_gtf(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    data = [line.strip().split('\t') for line in lines if 'FPKM' in line]

    result = []

    for line in data:
        transcript_id = re.search(r'transcript_id "([^"]*)"', line[8]).group(1)
        FPKM_value = re.search(r'FPKM "([^"]*)"', line[8]).group(1)
        result.append((transcript_id, FPKM_value))

    return result

def main(file_path, output_file):
    output_name = os.path.splitext(os.path.basename(output_file))[0]
    with open(output_file, 'w') as out_file:
        out_file.write(f"#transcript_id\t{output_name}\n")
        data = process_gtf(file_path)
        for row in data:
            out_file.write(f"{row[0]}\t{row[1]}\n")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script.py input_gtf output_file")
        sys.exit(1)
    
    input_gtf = sys.argv[1]
    output_file = sys.argv[2]

    main(input_gtf, output_file)
