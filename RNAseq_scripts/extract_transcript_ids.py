import sys

def extract_transcript_id(file_path):
    transcript_ids = set()

    with open(file_path, 'r') as file:
        for line in file:
            if 'transcript_id' in line:
                transcript_id = line.split('transcript_id "')[1].split('"')[0]
                transcript_ids.add(transcript_id)

    return transcript_ids

def main(input_file, output_file):
    transcript_ids = extract_transcript_id(input_file)

    with open(output_file, 'w') as out_file:
        for transcript_id in transcript_ids:
            out_file.write(f"{transcript_id}\n")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script.py input_gtf output_file")
        sys.exit(1)
    
    input_gtf = sys.argv[1]
    output_file = sys.argv[2]

    main(input_gtf, output_file)
