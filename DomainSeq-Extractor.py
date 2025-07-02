#!/usr/bin/env python3

import sys
from pathlib import Path
from typing import List, Tuple, Dict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def load_sequences(fasta_file: Path) -> Dict[str, SeqRecord]:
    """Load sequences from FASTA file into dictionary."""
    try:
        return SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    except FileNotFoundError:
        print(f"Error: FASTA file '{fasta_file}' not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading FASTA file: {e}", file=sys.stderr)
        sys.exit(1)

def parse_coordinates(coord_file: Path) -> List[Tuple[str, int, int]]:
    """Parse coordinate file and return list of (protein_id, start, end) tuples."""
    coords = []
    try:
        with open(coord_file, "r") as f:
            next(f)  # Skip header
            for line_num, line in enumerate(f, start=2):
                line = line.strip()
                if not line:
                    continue
                
                try:
                    parts = line.split("\t")
                    if len(parts) != 2:
                        print(f"Warning: Invalid format at line {line_num}, skipping.", file=sys.stderr)
                        continue
                    
                    protein_id, coord_range = parts
                    if "-" not in coord_range:
                        print(f"Warning: Invalid coordinate format '{coord_range}' at line {line_num}, skipping.", file=sys.stderr)
                        continue
                    
                    start_str, end_str = coord_range.split("-", 1)
                    start, end = int(start_str), int(end_str)
                    
                    if start <= 0 or end <= 0 or start > end:
                        print(f"Warning: Invalid coordinates {start}-{end} at line {line_num}, skipping.", file=sys.stderr)
                        continue
                    
                    coords.append((protein_id.strip(), start, end))
                    
                except ValueError as e:
                    print(f"Warning: Error parsing line {line_num}: {e}, skipping.", file=sys.stderr)
                    continue
                    
    except FileNotFoundError:
        print(f"Error: Coordinate file '{coord_file}' not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading coordinate file: {e}", file=sys.stderr)
        sys.exit(1)
    
    return coords

def extract_domains(seq_dict: Dict[str, SeqRecord], 
                   coords: List[Tuple[str, int, int]], 
                   output_file: Path) -> None:
    """Extract domain sequences and write to output file."""
    successful_extractions = 0
    failed_extractions = 0
    
    try:
        with open(output_file, "w") as out:
            for protein_id, start, end in coords:
                if protein_id not in seq_dict:
                    print(f"Warning: Protein '{protein_id}' not found in FASTA file.", file=sys.stderr)
                    failed_extractions += 1
                    continue

                full_seq = seq_dict[protein_id].seq
                seq_length = len(full_seq)
                
                # Validate coordinates against sequence length
                if end > seq_length:
                    print(f"Warning: End coordinate {end} exceeds sequence length {seq_length} for {protein_id}, skipping.", file=sys.stderr)
                    failed_extractions += 1
                    continue
                
                # Extract domain (convert from 1-based to 0-based indexing)
                domain_seq = full_seq[start - 1:end]
                
                # Write to output file
                out.write(f">{protein_id}_{start}-{end}\n")
                out.write(f"{domain_seq}\n")
                successful_extractions += 1
                
    except Exception as e:
        print(f"Error writing to output file: {e}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Successfully extracted {successful_extractions} domain sequences.")
    if failed_extractions > 0:
        print(f"Failed to extract {failed_extractions} domain sequences (see warnings above).")

def main():
    """Main function to orchestrate domain extraction."""
    # Configuration
    fasta_file = Path("D:/ST_Protein Motif/STP.fasta")
    coord_file = Path("D:/ST_Protein Motif/STP_CORD.txt")
    output_file = Path("D:/ST_Protein Motif/STP_Domain.fasta")
    
    # Validate input files exist
    if not fasta_file.exists():
        print(f"Error: FASTA file '{fasta_file}' does not exist.", file=sys.stderr)
        sys.exit(1)
    
    if not coord_file.exists():
        print(f"Error: Coordinate file '{coord_file}' does not exist.", file=sys.stderr)
        sys.exit(1)
    
    # Create output directory if it doesn't exist
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    print("Loading sequences...")
    seq_dict = load_sequences(fasta_file)
    print(f"Loaded {len(seq_dict)} sequences.")
    
    print("Parsing coordinates...")
    coords = parse_coordinates(coord_file)
    print(f"Parsed {len(coords)} coordinate entries.")
    
    if not coords:
        print("No valid coordinates found. Exiting.", file=sys.stderr)
        sys.exit(1)
    
    print("Extracting domains...")
    extract_domains(seq_dict, coords, output_file)
    
    print(f"Domain sequences written to '{output_file}'")

if __name__ == "__main__":
    main()