# Input data
input_file = "miat.collinearity"  # Replace with your input file name
output_file = "gene_pairs.txt"  # Replace with your desired output file name

# Read and process the data
with open(input_file, "r") as file:
    lines = file.readlines()

gene_pairs = []
for line in lines:
    # Ignore non-data lines
    if line.strip().startswith("#") or line.strip() == "" or line.strip().startswith("##"):
        continue

    # Extract relevant gene pair lines
    if "\t" in line:
        parts = line.strip().split("\t")
        if len(parts) >= 3:  # Ensure the line has at least two columns
            gene1 = parts[1]
            gene2 = parts[2]
            gene_pairs.append(f"{gene1}\t{gene2}")

# Write the gene pairs to the output file
with open(output_file, "w") as file:
    file.write("\n".join(gene_pairs))

print(f"Gene pairs have been extracted and saved to {output_file}")
