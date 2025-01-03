# Open and read the data
input_file = "Mango_Sugar_Transproter/mi.gff"  # Replace with your input file name
output_file = "Mango_Sugar_Transproter/mi_TEST.gff"  # Replace with your desired output file name

# Read the data
with open(input_file, "r") as file:
    lines = file.readlines()

# Process the data
modified_lines = []
for line in lines:
    columns = line.strip().split("\t")
    columns[0] = "Michr" + columns[0]  # Add "Mi" to the first column
    modified_lines.append("\t".join(columns))

# Write the modified data to a new file
with open(output_file, "w") as file:
    file.write("\n".join(modified_lines))

print(f"Data added to the first column has been saved to {output_file}")
