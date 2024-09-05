import sys
import os
import pandas as pd

def merge_outputs(output_files):
    dfs = []

    for output_file in output_files:
        df = pd.read_csv(output_file, sep='\t')
        df.set_index('#transcript_id', inplace=True)  # Adjust column name here
        dfs.append(df)

    merged_df = pd.concat(dfs, axis=1)
    return merged_df

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python merge_outputs.py output1.tsv output2.tsv ... outputN.tsv")
        sys.exit(1)
    
    output_files = sys.argv[1:]
    
    # Merge the outputs
    merged_df = merge_outputs(output_files)

    # Save the merged DataFrame to a TSV file
    merged_df.to_csv('merged_output.txt', sep='\t')

    print("Merged output saved as 'merged_output.txt'.")
