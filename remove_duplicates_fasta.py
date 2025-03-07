from Bio import SeqIO

def remove_duplicates_fasta_by_header(input_file, output_file):
    headers = set()  # Set to store unique headers
    unique_records = []  # List to store unique records
    
    # Read the input FASTA file
    with open(input_file, "r") as infile:
        for record in SeqIO.parse(infile, "fasta"):
            header = record.id  # Get the header (sequence ID)
            if header not in headers:
                headers.add(header)  # Add header to the set
                unique_records.append(record)  # Store the unique record
    
    # Write unique records to the output FASTA file
    with open(output_file, "w") as outfile:
        SeqIO.write(unique_records, outfile, "fasta")
    
    print(f"Total unique headers: {len(headers)}")
    print(f"Output written to: {output_file}")

# Example usage
if __name__ == "__main__":
    input_file = "input.fasta"
    output_file = "output_unique_by_header.fasta"
    remove_duplicates_fasta_by_header(input_file, output_file)
