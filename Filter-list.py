import pandas as pd

# File paths
gff_file = "G:\Bioinformatics_Data\mango_sugar_transproter\Micro Synteny\Circos\gene_dicot.txt"  # Replace with your GFF file path
gene_ids_file = "G:\Bioinformatics_Data\mango_sugar_transproter\Micro Synteny\Circos\gene_highlight.txt"  # Replace with your gene IDs file path
output_file = "G:\Bioinformatics_Data\mango_sugar_transproter\Micro Synteny\Circos\highlight.csv"  # Output CSV file path

# Load the GFF file
gff_columns = ["Chr", "GeneID", "Start", "End"]
gff_df = pd.read_csv(gff_file, sep="\t", header=None, names=gff_columns)

# Load the gene IDs
with open(gene_ids_file, "r") as file:
    gene_ids = [line.strip() for line in file]

# Filter the GFF data for matching gene IDs
filtered_df = gff_df[gff_df["GeneID"].isin(gene_ids)]

# Save the filtered data to a CSV file
filtered_df.to_csv(output_file, index=False)

print(f"Filtered data saved to {output_file}")
