import csv

# Input file paths
gff_file = "G:\Bioinformatics_Data\mango_sugar_transproter\Micro Synteny\Circos\gene_dicot.txt"  # Replace with the actual GFF file name
gene_pair_file = "G:\Bioinformatics_Data\mango_sugar_transproter\Micro Synteny\Circos\gene_pairs_dicot.txt"  # Replace with the actual gene pair file name
output_file = "G:\Bioinformatics_Data\mango_sugar_transproter\Micro Synteny\Circos\output.csv"  # Output CSV file

# Step 1: Parse the GFF file and create a mapping of feature IDs to chromosome, start, and end positions
gff_data = {}
with open(gff_file, "r") as gff:
    for line in gff:
        parts = line.strip().split("\t")
        if len(parts) == 4:  # Ensure proper format
            chromosome, feature_id, start, end = parts[0], parts[1], int(parts[2]), int(parts[3])
            gff_data[feature_id] = (chromosome, start, end)

# Step 2: Process the gene pair file and fetch chromosome, start, and end positions for each gene
output_data = []
with open(gene_pair_file, "r") as pairs:
    for line in pairs:
        gene1, gene2 = line.strip().split("\t")
        gene1_data = gff_data.get(gene1, ("NA", "NA", "NA"))
        gene2_data = gff_data.get(gene2, ("NA", "NA", "NA"))
        output_data.append([
            gene1, gene1_data[0], gene1_data[1], gene1_data[2],  # Gene1: ID, Chromosome, Start, End
            gene2, gene2_data[0], gene2_data[1], gene2_data[2]   # Gene2: ID, Chromosome, Start, End
        ])

# Step 3: Write the results to a CSV file
with open(output_file, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    # Write header
    writer.writerow([
        "Gene1", "Gene1_Chromosome", "Gene1_Start", "Gene1_End",
        "Gene2", "Gene2_Chromosome", "Gene2_Start", "Gene2_End"
    ])
    # Write rows
    writer.writerows(output_data)

print(f"Output written to {output_file}")
