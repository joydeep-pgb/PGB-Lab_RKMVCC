def extract_gene_ids_from_fasta(fasta_file):
    gene_ids = []

    with open(fasta_file, 'r') as file:
        lines = file.readlines()

        for line in lines:
            if line.startswith('>'):  # Identifying the header line
                gene_id = line.strip().split()[0][1:]  # Extracting the gene ID
                gene_ids.append(gene_id)

    return gene_ids

# Replace 'your_fasta_file.fasta' with the path to your actual FASTA file
fasta_file_path = 'DS_CDS.fasta'

# Extracting gene IDs
gene_ids_list = extract_gene_ids_from_fasta(fasta_file_path)

# Specify the path for the output text file
output_file_path = 'DS_DEG.txt'

# Writing gene IDs to a text file
with open(output_file_path, 'w') as output_file:
    for gene_id in gene_ids_list:
        output_file.write(f'{gene_id}\n')

print(f'Gene IDs extracted from {fasta_file_path} and saved to {output_file_path}')
