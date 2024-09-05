def calculate_gc_content(sequence):
    gc_count = sequence.count('G') + sequence.count('C')
    total_count = len(sequence)
    gc_content = (gc_count / total_count) * 100
    return gc_content

def read_fasta(file_path):
    sequences = {}
    with open(file_path, 'r') as file:
        lines = file.readlines()
        current_sequence = ''
        for line in lines:
            if line.startswith('>'):
                if current_sequence:
                    sequences[header] = current_sequence
                header = line.strip()[1:]
                current_sequence = ''
            else:
                current_sequence += line.strip()
        if current_sequence:
            sequences[header] = current_sequence
    return sequences

def main(file_path, output_file):
    sequences = read_fasta(file_path)
    with open(output_file, 'w') as outfile:
        outfile.write("Header\tGC Content (%)\n")
        for header, sequence in sequences.items():
            gc_content = calculate_gc_content(sequence)
            outfile.write(f"{header}\t{gc_content:.2f}\n")

if __name__ == "__main__":
    file_path = "/media/pgblab/Elements/RNAseq_Analysis/mRNA_seq.fasta"  # Replace with the path to your FASTA file
    output_file = "gc_content_results.tsv"  # Define the output file name
    main(file_path, output_file)

