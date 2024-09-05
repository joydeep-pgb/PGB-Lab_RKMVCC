def read_fasta_file(file_path):
    sequences = {}
    with open(file_path, 'r') as file:
        lines = file.readlines()
        current_sequence_id = ''
        for line in lines:
            if line.startswith('>'):
                current_sequence_id = line.strip()[1:]
                sequences[current_sequence_id] = ''
            else:
                sequences[current_sequence_id] += line.strip()
    return sequences

def extract_genes_from_fasta(gene_list_file, column_index, fasta_file):
    # Read gene IDs from the gene list file and extract specified column
    gene_ids = set()
    with open(gene_list_file, 'r') as gene_file:
        lines = gene_file.readlines()
        for line in lines:
            columns = line.strip().split()  # Split the line into columns
            if len(columns) > column_index:  # Check if the specified column index exists
                gene_ids.add(columns[column_index])  # Add gene ID to set

    # Read sequences from the FASTA file
    fasta_sequences = read_fasta_file(fasta_file)

    # Extract sequences for the unique gene IDs
    extracted_sequences = {}
    for gene_id in gene_ids:
        if gene_id in fasta_sequences:
            extracted_sequences[gene_id] = fasta_sequences[gene_id]

    return extracted_sequences

def write_extracted_sequences_to_file(extracted_sequences, output_file):
    with open(output_file, 'w') as out_file:
        for gene_id, sequence in extracted_sequences.items():
            out_file.write(f'>{gene_id}\n{sequence}\n')

# Replace these file paths with your actual file paths
gene_list_file_path = '12h_DELs.txt'
fasta_file_path = '12h_lncRNA.fasta'
output_file_path = '12h_DELs.fasta'
column_index_to_use = 0  # Specify the column index (0-based) containing gene IDs

# Extract genes from the FASTA file based on the gene list and write to output file
extracted_genes = extract_genes_from_fasta(gene_list_file_path, column_index_to_use, fasta_file_path)
write_extracted_sequences_to_file(extracted_genes, output_file_path)
