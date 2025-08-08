import pandas as pd
import os
from pathlib import Path
import numpy as np
from functools import reduce
import gc
import sys

def optimize_dtypes(df, common_column):
    """Optimize data types to reduce memory usage"""
    for col in df.columns:
        if col == common_column:
            continue  # Keep common column as string for merging
        
        # Try to convert to numeric if possible
        try:
            # Check if all values are numeric (including "0" from fillna)
            numeric_values = pd.to_numeric(df[col], errors='coerce')
            if not numeric_values.isna().all():
                # If mostly integers, use int32
                if numeric_values.dropna().apply(lambda x: x.is_integer()).all():
                    df[col] = numeric_values.astype('Int32')
                else:
                    df[col] = numeric_values.astype('float32')
        except:
            # Keep as string/object if conversion fails
            pass
    
    return df

def merge_tables_simple(folder_path, common_column=None, output_file=None):
    """
    Simple and efficient merge of tables without external dependencies
    
    Parameters:
    folder_path (str): Path to folder containing .txt files
    common_column (str): Name of common column to merge on
    output_file (str): Output file path
    """
    
    folder_path = Path(folder_path)
    txt_files = list(folder_path.glob("*.txt"))
    
    if not txt_files:
        raise ValueError("No .txt files found")
    
    print(f"Found {len(txt_files)} files to merge")
    
    # Detect file format from first file
    first_file = txt_files[0]
    print(f"Analyzing file format using: {first_file.name}")
    
    # Try different separators
    separators = ['\t', ',', '|', ';']
    best_sep = '\t'
    max_cols = 0
    
    for sep in separators:
        try:
            sample_df = pd.read_csv(first_file, sep=sep, nrows=5, dtype=str)
            if len(sample_df.columns) > max_cols:
                max_cols = len(sample_df.columns)
                best_sep = sep
        except:
            continue
    
    print(f"Using separator: '{best_sep}' ({max_cols} columns detected)")
    
    # Read first file to understand structure
    first_df = pd.read_csv(first_file, sep=best_sep, dtype=str)
    first_df.columns = first_df.columns.str.strip()
    
    if common_column is None:
        # Auto-detect common column by reading a few more files
        common_columns = set(first_df.columns)
        
        # Check first 3 files to find truly common columns
        for file_path in txt_files[1:4]:
            try:
                temp_df = pd.read_csv(file_path, sep=best_sep, nrows=1, dtype=str)
                temp_df.columns = temp_df.columns.str.strip()
                common_columns &= set(temp_df.columns)
            except:
                continue
        
        if not common_columns:
            raise ValueError("No common columns found across files")
        
        common_column = list(common_columns)[0]
    
    print(f"Using common column: '{common_column}'")
    print(f"Columns in dataset: {list(first_df.columns)}")
    
    # Verify common column exists
    if common_column not in first_df.columns:
        raise ValueError(f"Column '{common_column}' not found in files")
    
    # Start with first file
    result_df = first_df.copy()
    print(f"Started with {first_file.name}: {len(result_df):,} rows")
    
    # Process remaining files
    for i, file_path in enumerate(txt_files[1:], 1):
        print(f"Processing file {i+1}/{len(txt_files)}: {file_path.name}")
        
        try:
            # Read file
            df = pd.read_csv(file_path, sep=best_sep, dtype=str)
            df.columns = df.columns.str.strip()
            
            print(f"  Loaded: {len(df):,} rows, {len(df.columns)} columns")
            
            # Verify common column exists
            if common_column not in df.columns:
                print(f"  Warning: Column '{common_column}' not found in {file_path.name}, skipping...")
                continue
            
            # Merge with result
            print(f"  Merging...")
            result_df = pd.merge(result_df, df, on=common_column, how='outer', suffixes=('', '_dup'))
            
            # Remove duplicate columns
            dup_cols = [col for col in result_df.columns if col.endswith('_dup')]
            if dup_cols:
                result_df = result_df.drop(columns=dup_cols)
            
            print(f"  Result: {len(result_df):,} rows, {len(result_df.columns)} columns")
            
            # Clean up
            del df
            gc.collect()
            
        except Exception as e:
            print(f"  Error processing {file_path.name}: {e}")
            continue
    
    print(f"\nMerge completed: {len(result_df):,} rows, {len(result_df.columns)} columns")
    
    # Fill empty cells with "0"
    print("Filling empty cells with '0'...")
    result_df = result_df.fillna("0")
    
    # Optimize data types to save memory
    print("Optimizing data types...")
    result_df = optimize_dtypes(result_df, common_column)
    
    # Display column info
    print(f"\nFinal dataset:")
    print(f"  Rows: {len(result_df):,}")
    print(f"  Columns: {len(result_df.columns)}")
    print(f"  Column names: {list(result_df.columns)}")
    
    # Show memory usage per column
    print(f"\nMemory usage by column:")
    memory_usage = result_df.memory_usage(deep=True)
    for col in result_df.columns:
        mb_usage = memory_usage[col] / (1024 * 1024)
        print(f"  {col}: {mb_usage:.2f} MB")
    
    total_mb = memory_usage.sum() / (1024 * 1024)
    print(f"  Total: {total_mb:.2f} MB")
    
    # Save result
    if output_file:
        print(f"\nSaving to {output_file}...")
        
        try:
            if output_file.endswith('.csv'):
                result_df.to_csv(output_file, index=False)
            elif output_file.endswith(('.xlsx', '.xls')):
                result_df.to_excel(output_file, index=False)
            else:
                # Default to CSV
                result_df.to_csv(output_file + '.csv', index=False)
                output_file = output_file + '.csv'
            
            print(f"Successfully saved to {output_file}")
            
            # Verify saved file
            saved_size_mb = Path(output_file).stat().st_size / (1024 * 1024)
            print(f"File size: {saved_size_mb:.2f} MB")
            
        except Exception as e:
            print(f"Error saving file: {e}")
    
    return result_df

def quick_preview(df, n_rows=5):
    """Show a preview of the merged data"""
    print(f"\nPreview of merged data (first {n_rows} rows):")
    print("=" * 80)
    
    # Show column names
    print("Columns:", list(df.columns))
    print()
    
    # Show first few rows
    preview = df.head(n_rows)
    
    # Format output to fit screen better
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)
    pd.set_option('display.max_colwidth', 15)
    
    print(preview.to_string(index=False))
    
    # Reset display options
    pd.reset_option('display.max_columns')
    pd.reset_option('display.width')  
    pd.reset_option('display.max_colwidth')

def main():
    """Main function"""
    
    # Configuration - update these paths
    folder_path = "F:\\ABC Transporters\\Phaseolus_NEW_DEG\\NEW_SRA\\Drought_Leaf\\Assembled\\TPM"  # Change this to your actual folder path
    common_column = "gene_id"  # Based on your previous error message
    output_file = "F:\\ABC Transporters\\Phaseolus_NEW_DEG\\NEW_SRA\\Drought_Leaf\\Assembled\\TPM\\merged_gene_data.csv"
    
    print("Starting table merge process...")
    print("=" * 50)
    
    try:
        # Perform the merge
        merged_df = merge_tables_simple(
            folder_path=folder_path,
            common_column=common_column,
            output_file=output_file
        )
        
        # Show preview
        quick_preview(merged_df)
        
        print("\n" + "="*50)
        print("MERGE COMPLETED SUCCESSFULLY!")
        print("="*50)
        print(f"Output file: {output_file}")
        print(f"Final dimensions: {merged_df.shape[0]:,} rows × {merged_df.shape[1]} columns")
        
        return merged_df
        
    except Exception as e:
        print(f"\nError during merge operation: {e}")
        print("\nDebugging suggestions:")
        print("1. Check that all .txt files have the same separator")
        print("2. Verify the common column name exists in all files")
        print("3. Check file paths are correct")
        print("4. Look for any corrupted or differently formatted files")
        
        return None

# Additional utility function for troubleshooting
def debug_files(folder_path, max_files=5):
    """Debug file structure to identify issues"""
    folder_path = Path(folder_path)
    txt_files = list(folder_path.glob("*.txt"))
    
    print("DEBUGGING FILE STRUCTURE")
    print("=" * 50)
    
    for i, file_path in enumerate(txt_files[:max_files]):
        print(f"\nFile {i+1}: {file_path.name}")
        print("-" * 40)
        
        try:
            # Try different separators
            for sep in ['\t', ',', '|', ';']:
                try:
                    df = pd.read_csv(file_path, sep=sep, nrows=3, dtype=str)
                    if len(df.columns) > 1:
                        print(f"Separator '{sep}': {len(df.columns)} columns")
                        df.columns = df.columns.str.strip()
                        print(f"Columns: {list(df.columns)}")
                        print("Sample data:")
                        print(df.head(2).to_string(index=False))
                        break
                except:
                    continue
                    
        except Exception as e:
            print(f"Error reading file: {e}")

if __name__ == "__main__":
    # Uncomment the line below to debug file structure if needed
    # debug_files("path/to/your/txt/files")
    
    main()