import pandas as pd
import os

def merge_tsv_files(folder_path, ref_column):
    # List all files in the given folder
    tsv_files = [f for f in os.listdir(folder_path) if f.endswith('.txt')]
    
    # Initialize an empty DataFrame
    merged_df = pd.DataFrame()
    
    # Iterate over the list of files and merge them
    for file in tsv_files:
        file_path = os.path.join(folder_path, file)
        df = pd.read_csv(file_path, sep='\t')
        
        if merged_df.empty:
            merged_df = df
        else:
            merged_df = pd.merge(merged_df, df, on=ref_column, how='outer')
    
    # Fill missing values with "0"
    merged_df.fillna(0, inplace=True)
    
    return merged_df

# Set the folder path and reference column name
folder_path = r"F:\Sorghum_MetaDEG\Batch\count drought"
ref_column = 'gene_id'

# Merge the files
merged_df = merge_tsv_files(folder_path, ref_column)

# Save the merged dataframe to a new TSV file
merged_df.to_csv('Sorghum_Drought_Counts.txt', sep='\t', index=False)
