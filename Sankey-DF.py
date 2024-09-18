import pandas as pd

# Load the input data from the text file
file_path = 'sankey.txt'  # Replace with your actual file path
df = pd.read_csv(file_path, sep='\t')  # Assuming your file is tab-separated

# Create the transformation steps without dropping duplicates

# Step 1: Create the lncRNA -> miRNA mapping
lncRNA_miRNA = df[['lncRNA', 'miRNA']].copy()  # Use .copy() to avoid SettingWithCopyWarning
lncRNA_miRNA.loc[:, 'output'] = lncRNA_miRNA['lncRNA'] + ' [1] ' + lncRNA_miRNA['miRNA']

# Step 2: Create the miRNA -> Target mapping
miRNA_Target = df[['miRNA', 'Target']].copy()  # Use .copy() to avoid SettingWithCopyWarning
miRNA_Target.loc[:, 'output'] = miRNA_Target['miRNA'] + ' [1] ' + miRNA_Target['Target']

# Step 3: Create the Target -> GO mapping
Target_GO = df[['Target', 'GO']].copy()  # Use .copy() to avoid SettingWithCopyWarning
Target_GO.loc[:, 'output'] = Target_GO['Target'] + ' [1] ' + Target_GO['GO']

# Combine all outputs into a single list
combined_output = pd.concat([lncRNA_miRNA['output'], miRNA_Target['output'], Target_GO['output']])

# Save the final output to a new text file
output_file = 'sankey_output.txt'  # Specify your desired output file path
with open(output_file, 'w') as f:
    for line in combined_output:
        f.write(line + '\n')

print(f'Output saved to {output_file}')
