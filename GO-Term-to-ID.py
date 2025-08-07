import pandas as pd

# Sample data - replace this with your actual data loading method
data = pd.read_csv("F:\\Mango_AS_rMATS\\rMATS\\DAS_GO-KEGG\\GFvsRF_GO.txt", sep='\t')

# Split the 'Objects' column into a list and explode it
data['Objects'] = data['Objects'].apply(lambda x: x.split(','))

# Explode the DataFrame to create a row for each object
data_exploded = data.explode('Objects')

# Display the modified DataFrame
print(data_exploded)

# Save the modified DataFrame to a text file
data_exploded.to_csv("F:\\Mango_AS_rMATS\\rMATS\\DAS_GO-KEGG\\GFvsRF_GO.txt", sep='\t', index=False)

# Optional: Print confirmation
print("Data has been Processed and saved successfully.")
