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
        
        # Get filename without extension
        filename_without_ext = os.path.splitext(file)[0]
        
        # Keep only Gene ID, TPM, and FPKM columns
        columns_to_keep = [ref_column]
        if 'TPM' in df.columns:
            columns_to_keep.append('TPM')
        if 'FPKM' in df.columns:
            columns_to_keep.append('FPKM')
        
        df_filtered = df[columns_to_keep]
        
        if merged_df.empty:
            merged_df = df_filtered
        else:
            # Create a copy of filtered df
            df_copy = df_filtered.copy()
            
            # Rename TPM and FPKM columns with filename prefix
            new_columns = {}
            for col in df_copy.columns:
                if col != ref_column:
                    new_columns[col] = f"{filename_without_ext}_{col}"
            
            df_copy.rename(columns=new_columns, inplace=True)
            
            # Merge with the main dataframe
            merged_df = pd.merge(merged_df, df_copy, on=ref_column, how='outer')
    
    # Fill missing values with "0"
    merged_df.fillna(0, inplace=True)
    
    return merged_df

# Set the folder path and reference column name
folder_path = "G:\Sorghum_MetaDEG\Drought\Assembled\TPMs"
ref_column = "Gene ID"

# Merge the files
merged_df = merge_tsv_files(folder_path, ref_column)

# Save the merged dataframe to a new TSV file
merged_df.to_csv('G:\Sorghum_MetaDEG\Drought\Assembled\TPMs\Sbi_DR_Counts.txt', sep='\t', index=False)