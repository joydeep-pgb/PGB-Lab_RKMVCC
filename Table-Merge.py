import pandas as pd

# Function to merge two DataFrames based on a single column with a specified join type
def merge_dataframes(file1, file2, output_file, join_column='ID', join_type='inner'):
    """
    Merge two DataFrames from text files based on a specified column and save the result to a text file.
    
    Parameters:
    - file1: Path to the first input text file.
    - file2: Path to the second input text file.
    - output_file: Path to the output text file.
    - join_column: Column to join on (default is 'ID').
    - join_type: Type of join ('inner', 'outer', 'left', 'right') (default is 'inner').
    
    Join Types:
    - 'inner': Returns only the rows with matching keys in both DataFrames.
    - 'outer': Returns all rows from both DataFrames, with NaNs in places where one DataFrame is missing data.
    - 'left': Returns all rows from the left DataFrame and the matching rows from the right DataFrame. Non-matching rows from the right DataFrame result in NaNs.
    - 'right': Returns all rows from the right DataFrame and the matching rows from the left DataFrame. Non-matching rows from the left DataFrame result in NaNs.
    
    Returns:
    - None
    """
    # Load the DataFrames from the text files
    df1 = pd.read_csv(file1, sep='\t')
    df2 = pd.read_csv(file2, sep='\t')
    
    # Merge the DataFrames on the specified column
    merged_df = pd.merge(df1, df2, on=join_column, how=join_type)
    
    # Save the merged DataFrame to a text file
    merged_df.to_csv(output_file, sep='\t', index=False)
    
    # Optional: Print confirmation
    print(f"Merged data has been saved to '{output_file}'.")

# Example usage of the function
merge_dataframes('Table1.txt', 'Table2.txt', 'Little-fruit_STATs.txt', join_column='Sample Name', join_type='inner')
