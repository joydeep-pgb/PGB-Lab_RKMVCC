import csv

# Input file paths
gff_file = "G:\Bioinformatics_Data\mango_sugar_transproter\Micro Synteny\Circos\gene_dicot.txt"  # Replace with the actual GFF file name
gene_pair_file = "G:\Bioinformatics_Data\mango_sugar_transproter\Micro Synteny\Circos\gene_pairs_dicot.txt"  # Replace with the actual gene pair file name
output_file = "G:\Bioinformatics_Data\mango_sugar_transproter\Micro Synteny\Circos\output.txt"  # Output CSV file

# Step 1: Parse the GFF file and create a mapping of feature IDs to start and end positions
gff_data = {}
with open(gff_file, "r") as gff:
    for line in gff:
        parts = line.strip().split("\t")
        if len(parts) == 4:  # Ensure proper format
            feature_id, start, end = parts[1], int(parts[2]), int(parts[3])
            gff_data[feature_id] = (start, end)

# Step 2: Process the gene pair file and fetch positions for each gene
output_data = []
with open(gene_pair_file, "r") as pairs:
    for line in pairs:
        gene1, gene2 = line.strip().split("\t")
        gene1_start, gene1_end = gff_data.get(gene1, ("NA", "NA"))
        gene2_start, gene2_end = gff_data.get(gene2, ("NA", "NA"))
        output_data.append([gene1, gene1_start, gene1_end, gene2, gene2_start, gene2_end])

# Step 3: Write the results to a CSV file
with open(output_file, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["Gene1", "Gene1_Start", "Gene1_End", "Gene2", "Gene2_Start", "Gene2_End"])  # Header
    writer.writerows(output_data)

print(f"Output written to {output_file}")
