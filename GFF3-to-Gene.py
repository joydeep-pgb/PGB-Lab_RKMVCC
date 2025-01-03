import pandas as pd

# Input GFF file
input_file = "Osativa_323_v7.0.gene.gff3"
output_file = "OS.gtf"

# Columns to extract: Chromosome (1st), Gene Name (Name in the 9th column), Start (4th), End (5th)
data = []
with open(input_file, 'r') as file:
    for line in file:
        if line.startswith('#'):
            continue  # Skip header lines
        columns = line.strip().split('\t')
        if columns[2] == 'gene':  # Only extract "gene" rows
            chr = columns[0]
            start = columns[3]
            end = columns[4]
            attributes = columns[8]
            # Extract gene name from attributes
            gene_name = [attr.split('=')[1] for attr in attributes.split(';') if attr.startswith('Name=')][0]
            data.append([chr, gene_name, start, end])

# Convert to DataFrame and save as TSV
df = pd.DataFrame(data, columns=['Chromosome', 'Gene_Name', 'Start', 'End'])
df.to_csv(output_file, sep='\t', index=False)

print(f"Extracted data saved to {output_file}")
