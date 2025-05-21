from Bio import SeqIO

def read_fasta_file(file_path):
    """Reads a FASTA file into a dictionary."""
    return {record.id: str(record.seq) for record in SeqIO.parse(file_path, "fasta")}

def extract_genes_from_fasta(gene_ids, fasta_file):
    """Extracts sequences matching given gene IDs."""
    fasta_sequences = read_fasta_file(fasta_file)
    return {gene_id: fasta_sequences[gene_id] for gene_id in gene_ids if gene_id in fasta_sequences}

def write_extracted_sequences_to_file(extracted_sequences, output_file, line_length=60):
    """Writes sequences to a FASTA file with line wrapping."""
    with open(output_file, 'w') as out_file:
        for gene_id, sequence in extracted_sequences.items():
            out_file.write(f'>{gene_id}\n')
            for i in range(0, len(sequence), line_length):
                out_file.write(sequence[i:i+line_length] + '\n')

# Paste your plain text IDs here (keep the triple quotes)
plain_text_ids = """
Manin16g008780.1
Manin00g005170.1
Manin07g001560.1
Manin08g008090.1
Manin04g012780.1
Manin20g006700.1
Manin09g014240.1
Manin02g003920.1
Manin16g003480.1
Manin09g013320.1
Manin18g002090.1
Manin03g006160.1
Manin01g010470.1
Manin10g005510.1
Manin17g008160.1
Manin00g002200.1
Manin00g002780.1
Manin12g010140.1
Manin03g002300.1
Manin07g003980.1
Manin19g004240.1
"""

# Convert plain text to Python string set
gene_ids = {line.strip() for line in plain_text_ids.splitlines() if line.strip()}

# File paths (update these as needed)
fasta_file_path = r"D:\Mango Genome_Annotation\mango_TA4_proteins_renamed.fasta"
output_file_path = r"C:\Users\ADMIN\Desktop\Documents\Mango ABC\FC1\UFvsRF_ABC.fasta"

# Execute the workflow
fasta_sequences = read_fasta_file(fasta_file_path)
extracted_genes = extract_genes_from_fasta(gene_ids, fasta_file_path)
write_extracted_sequences_to_file(extracted_genes, output_file_path)

# Print results
print(f"Total IDs requested: {len(gene_ids)}")
print(f"Successfully extracted: {len(extracted_genes)}")
if len(gene_ids) != len(extracted_genes):
    missing = gene_ids - extracted_genes.keys()
    print(f"Missing IDs ({len(missing)}):")
    for id in sorted(missing):
        print(f" - {id}")