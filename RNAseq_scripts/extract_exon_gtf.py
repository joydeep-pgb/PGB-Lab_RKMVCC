import re

def extract_gtf_info(gtf_line):
    # Define regular expressions for extracting information
    fields = gtf_line.strip().split('\t')

    if len(fields) < 9:
        return None

    attributes = fields[8]
    
    # Extract information using string splitting
    transcript_id_match = re.search(r'transcript_id "([^"]+)"', attributes)
    exon_number_match = re.search(r'exon_number "([^"]+)"', attributes)

    # Check if all necessary information is found
    if transcript_id_match and exon_number_match:
        chromosome = fields[0]
        transcript_id = transcript_id_match.group(1)
        start = int(fields[3])
        end = int(fields[4])
        exon_number = exon_number_match.group(1)

        return chromosome, transcript_id, start, end, exon_number

    return None

# Read the GTF file and write results to a text file in TSV format
gtf_file_path = '12h_sbi_mergeanno.gtf'  # Replace with the actual file path
output_file_path = 'output.tsv'  # Replace with the desired output file path

with open(gtf_file_path, 'r') as gtf_file, open(output_file_path, 'w') as output_file:
    # Write header to the output file
    output_file.write("Chromosome\tTranscript_ID\tStart\tEnd\tExon_Number\n")

    # Process each line in the GTF file
    for line in gtf_file:
        # Skip lines that do not contain 'transcript' or 'exon'
        if 'transcript' in line or 'exon' in line:
            result = extract_gtf_info(line)
            if result:
                chromosome, transcript_id, start, end, exon_number = result
                output_file.write(f"{chromosome}\t{transcript_id}\t{start}\t{end}\t{exon_number}\n")

print(f"Results written to {output_file_path}")
