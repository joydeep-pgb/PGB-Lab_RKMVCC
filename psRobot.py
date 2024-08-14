# Function to remove the first newline after the FASTA header and make it a single line while retaining tabs
def remove_first_newline_after_header(fasta_text):
    lines = fasta_text.splitlines()
    result = []
    for i, line in enumerate(lines):
        if line.startswith(">") and i < len(lines) - 1:
            result.append(line + lines[i + 1])  # Append the header and next line as a single line
        elif not line.startswith(">"):
            if not lines[i - 1].startswith(">"):
                result.append(line)  # Append other lines normally
        else:
            pass  # Skip the next line as it has already been appended
    return "\n".join(result)

# Read input from a text file
input_file = "MangoMeta_lncRNA/mango_miRNA-target.gTP"
with open(input_file, "r") as file:
    fasta_text = file.read()

# Process the FASTA text
processed_fasta_text = remove_first_newline_after_header(fasta_text)

# Write the output to a new text file
output_file = "MangoMeta_lncRNA/mango_miRNA-target.txt"
with open(output_file, "w") as file:
    file.write(processed_fasta_text)

print(f"Processed FASTA text has been saved to {output_file}")
