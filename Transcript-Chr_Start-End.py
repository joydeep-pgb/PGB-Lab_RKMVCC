# Open the GTF file for reading
input_file = 'MangoMeta_lncRNA/Mango_merge.gtf'  # Replace with your actual GTF file path
output_file = 'transcripts_info.txt'

# Initialize an empty list to store the extracted transcript information
transcripts = []

# Open and read the GTF file
with open(input_file, 'r') as gtf:
    for line in gtf:
        # Skip comments or headers
        if line.startswith("#"):
            continue

        # Split the GTF line into fields
        fields = line.strip().split('\t')

        # Check if the line describes a transcript (and not an exon)
        if fields[2] == "transcript":
            # Extract relevant fields
            chr_num = fields[0]          # Chromosome number
            start = fields[3]            # Transcript start position
            end = fields[4]              # Transcript end position
            # Extract transcript ID from the attributes field (9th column)
            attributes = fields[8]
            transcript_id = attributes.split('transcript_id "')[1].split('"')[0]

            # Add the extracted data to the transcripts list
            transcripts.append(f"{chr_num}\t{transcript_id}\t{start}\t{end}")

# Write the extracted data to the output file
with open(output_file, 'w') as out:
    out.write("Chr\tTranscript_ID\tStart\tEnd\n")  # Add a header
    for transcript in transcripts:
        out.write(transcript + "\n")

print(f"Transcript information has been saved to {output_file}")
