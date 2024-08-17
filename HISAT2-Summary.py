import os
import re
import pandas as pd

# Define the directory containing the text files
directory = 'D:\Bioinformatics_Data\Mango Meta analysis\Mapped Reads\Summary\Developing fruit Mapping summary'  # Replace with the path to your folder

# Initialize a list to hold the data for the DataFrame
data = []

# Regular expression pattern to match the "Aligned concordantly 1 time" line
pattern = re.compile(r'Aligned concordantly 1 time:\s*([\d.]+) \(([\d.]+)%\)')

# Iterate over each file in the directory
for filename in os.listdir(directory):
    if filename.endswith('.txt'):
        sample_name = os.path.splitext(filename)[0]
        file_path = os.path.join(directory, filename)

        with open(file_path, 'r') as file:
            content = file.read()
            match = pattern.search(content)
            if match:
                percentage = match.group(2)  # Extract the percentage
                data.append([sample_name, percentage])

# Create a DataFrame from the collected data
df = pd.DataFrame(data, columns=['Sample Name', 'Uniquely mapped (%)'])

# Define the output file
output_file = 'alignment_summary.txt'

# Write the DataFrame to a CSV file
df.to_csv(output_file, sep='\t', index=False)

print(f"Results have been written to {output_file}")
