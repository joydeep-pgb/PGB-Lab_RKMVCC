# Define the file name
input_file = "miRNA_data.txt"

# Read the content of the file
with open(input_file, 'r') as file:
    data = file.read()

# Split data by "//" to separate each pair
pairs = data.strip().split("//")

# Initialize a list to store pairs with gap = 3 and a list for miRNA-target IDs
filtered_pairs = []
miRNA_target_list = []

# Loop through each pair
for pair in pairs:
    # Check if the pair contains 'gap' and its value is 3
    if "gap       4" in pair:
        filtered_pairs.append(pair.strip())
        
        # Extract miRNA and target IDs
        lines = pair.strip().split("\n")
        miRNA = lines[0].split()[1]  # Extract miRNA ID
        target = lines[1].split()[1]  # Extract target ID
        miRNA_target_list.append((miRNA, target))

# Join the filtered pairs back with "\n//" as a separator to add a line break before "//"
result = "\n//\n".join(filtered_pairs) + "\n//"

# Print the filtered pairs
print("Filtered miRNA-target pairs:\n")
print(result)

# Print the miRNA and target IDs in tab-separated format
print("\nmiRNA and Target IDs (where gap = 3):\n")
print("miRNA ID\tTarget ID")
for miRNA, target in miRNA_target_list:
    print(f"{miRNA}\t{target}")

# Optionally, you can save the result to a new file
output_file = "filtered_miRNA_data.txt"
with open(output_file, 'w') as file:
    file.write(result)

# Optionally, save the tab-separated miRNA-target IDs to a file
output_table_file = "miRNA_target_pairs.txt"
with open(output_table_file, 'w') as file:
    file.write("miRNA ID\tTarget ID\n")
    for miRNA, target in miRNA_target_list:
        file.write(f"{miRNA}\t{target}\n")
