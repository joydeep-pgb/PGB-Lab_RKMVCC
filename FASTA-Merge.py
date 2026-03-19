from Bio import SeqIO
from pathlib import Path
import itertools

# 1. Set your main directory and the output file name
base_dir = Path(r"G:\Sorghum_circRNA\Salt\CIRI-full_Complete\Results\circRNA_Sequences")
output_file = base_dir / "Salt_circRNAs.fa"

def merge_fasta_files():
    all_sequences = []
    file_count = 0
    
    print("Searching for .fa and .fasta files...\n")
    
    # 2. Search for both extensions using wildcards
    fa_files = base_dir.rglob("*.fa")
    fasta_files = base_dir.rglob("*.fasta")
    
    # Combine the search results using itertools
    all_target_files = itertools.chain(fa_files, fasta_files)
    
    # 3. Process each found file
    for fasta_path in all_target_files:
        # Safety Check: Skip the output file itself so we don't create an infinite loop
        if fasta_path.name == output_file.name:
            continue
            
        file_count += 1
        
        # Read the sequences from the current file
        records = list(SeqIO.parse(fasta_path, "fasta"))
        all_sequences.extend(records)
        
        # Print progress to see which file was just processed
        print(f"Processed {len(records)} sequences from: {fasta_path.name}")

    # 4. Write everything into one single FASTA file
    if all_sequences:
        SeqIO.write(all_sequences, output_file, "fasta")
        print("\n" + "-"*40)
        print(f"Success! Merged {len(all_sequences)} total sequences from {file_count} files.")
        print(f"Saved to: {output_file}")
    else:
        print("No .fa or .fasta files were found. Please check your folder path.")

if __name__ == "__main__":
    merge_fasta_files()