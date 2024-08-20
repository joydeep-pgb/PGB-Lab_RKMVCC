import pandas as pd

def parse_input_file(input_file):
    # Initialize lists to hold the data
    gene_ids = []
    forward_primers = []
    reverse_primers = []
    amplicon_lengths = []

    # Read the input file line by line
    with open(input_file, 'r') as file:
        lines = file.readlines()[1:]  # Skip header line

        # Iterate through the lines to capture each set of data
        for i in range(0, len(lines), 2):
            # First line contains Gene ID, Forward primer, and Amplicon length
            line1 = lines[i].strip().split('\t')
            
            # Ensure the line is properly formatted
            if len(line1) == 3:
                gene_id = line1[0]
                forward_primer = line1[1].split('F: ')[1]
                amplicon_length = line1[2]
            else:
                print(f"Warning: Unexpected format in line {i+2}: {lines[i]}")
                continue
            
            # Second line contains only the Reverse primer
            reverse_primer = lines[i + 1].strip().split('R: ')[1]
            
            # Append the extracted data to the lists
            gene_ids.append(gene_id)
            forward_primers.append(forward_primer)
            reverse_primers.append(reverse_primer)
            amplicon_lengths.append(amplicon_length)

    # Create a DataFrame from the collected data
    data = {
        "Gene ID": gene_ids,
        "Forward primer sequence（5'→3'）": forward_primers,
        "Reverse primer sequence（5'→3'）": reverse_primers,
        "Amplicon length": amplicon_lengths
    }
    df = pd.DataFrame(data)
    return df

def save_to_file(df, output_file):
    df.to_csv(output_file, sep='\t', index=False)

def main():
    input_file = 'mango_primer.txt'  # Name of the input file
    output_file = 'Sample_output.txt'  # Name of the output file

    df = parse_input_file(input_file)
    if not df.empty:
        save_to_file(df, output_file)
        print(f"Data successfully saved to {output_file}")
    else:
        print("No valid data found to save.")

if __name__ == "__main__":
    main()
