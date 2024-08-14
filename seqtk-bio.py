from Bio import SeqIO

def read_fasta_file(file_path):
    sequences = {}
    for record in SeqIO.parse(file_path, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences

def extract_genes_from_fasta(gene_list_file, column_index, fasta_file):
    gene_ids = set()
    with open(gene_list_file, 'r') as gene_file:
        lines = gene_file.readlines()
        for line in lines:
            columns = line.strip().split()
            if len(columns) > column_index:
                gene_ids.add(columns[column_index])

    fasta_sequences = read_fasta_file(fasta_file)

    extracted_sequences = {}
    for gene_id in gene_ids:
        if gene_id in fasta_sequences:
            extracted_sequences[gene_id] = fasta_sequences[gene_id]

    return extracted_sequences

def write_extracted_sequences_to_file(extracted_sequences, output_file, line_length=60):
    with open(output_file, 'w') as out_file:
        for gene_id, sequence in extracted_sequences.items():
            out_file.write(f'>{gene_id}\n')
            for i in range(0, len(sequence), line_length):
                out_file.write(sequence[i:i+line_length] + '\n')

# Replace these file paths with your actual file paths
gene_list_file_path = 'MangoMeta_lncRNA/Order_ID.txt'
fasta_file_path = 'MangoMeta_lncRNA/mango_TA4_cds.fasta'
output_file_path = 'MangoMeta_lncRNA/qPCR_Target.fasta'
column_index_to_use = 0  # Specify the column index (0-based) containing gene IDs

# Extract genes from the FASTA file based on the gene list and write to output file
extracted_genes = extract_genes_from_fasta(gene_list_file_path, column_index_to_use, fasta_file_path)
write_extracted_sequences_to_file(extracted_genes, output_file_path)
