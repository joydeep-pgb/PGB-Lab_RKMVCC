import pandas as pd

# Load the exploded data from the text file
data_exploded = pd.read_csv("F:\\Mango_AS_rMATS\\rMATS\\DAS_GO-KEGG\\LFvsRF_GO.txt", sep='\t')

# Group by the relevant columns and join the 'Objects' values back together
data_original = data_exploded.groupby(['GO Term', 'Ontology', 'Description', 'Gene No', 'P-value'], as_index=False).agg({'Genes': ', '.join})

# Save the re-constructed DataFrame back to a text file
data_original.to_csv("F:\\Mango_AS_rMATS\\rMATS\\DAS_GO-KEGG\\LFvsRF_GO.txt", sep='\t', index=False)

# Optional: Print confirmation
print("Data has been reconstructed and saved.")
