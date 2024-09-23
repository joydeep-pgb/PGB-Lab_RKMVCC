import pandas as pd

# Load the exploded data from the text file
data_exploded = pd.read_csv('GO-Term.txt', sep='\t')

# Group by the relevant columns and join the 'genes' values back together
data_original = data_exploded.groupby(['GO Term', 'Ontology', 'Description', 'Gene No', 'pvalue'], as_index=False).agg({'genes': ','.join})

# Add a new column 'no gene' that counts the number of genes in the 'genes' column
data_original['no gene'] = data_original['genes'].apply(lambda x: len(x.split(',')))

# Save the re-constructed DataFrame back to a text file
data_original.to_csv('reconstructed_input.txt', sep='\t', index=False)

# Optional: Print confirmation
print("Data has been reconstructed, including the 'no gene' column, and saved to 'reconstructed_input.txt'.")
