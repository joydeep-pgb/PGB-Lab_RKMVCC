import sys

def extract_information(file_path):
    data = []

    with open(file_path, 'r') as file:
        for line in file:
            if 'transcript_id' in line:
                parts = line.strip().split('\t')
                chromosome = parts[0]
                start = parts[3]
                end = parts[4]
                exon_number = parts[8].split('exon_number "')[1].split('"')[0]
                transcript_id = parts[8].split('transcript_id "')[1].split('"')[0]
                data.append([chromosome, start, end, exon_number, transcript_id])

    return data

def save_to_tsv(data, output_file):
    with open(output_file, 'w') as file:
        file.write("Chromosome\tStart\tEnd\tExon_Number\tTranscript_ID\n")
        for row in data:
            file.write("\t".join(row) + "\n")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script.py input_gtf output_file")
        sys.exit(1)
    
    input_gtf = sys.argv[1]
    output_file = sys.argv[2]

    data = extract_information(input_gtf)
    save_to_tsv(data, output_file)

