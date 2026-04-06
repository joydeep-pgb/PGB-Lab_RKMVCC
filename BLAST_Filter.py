import sys

def filter_best_orthologs(input_file, output_file):
    # Dictionary to keep track of the best hit for each query gene
    best_hits = {}

    print(f"Reading from {input_file}...")
    
    with open(input_file, 'r') as infile:
        for line in infile:
            # Skip empty lines
            if not line.strip():
                continue
                
            parts = line.strip().split('\t')
            
            # Ensure the line has the expected 12 columns
            if len(parts) < 12:
                continue 

            query_id = parts[0]
            
            try:
                evalue = float(parts[10])
                bitscore = float(parts[11])
            except ValueError:
                # Skips any potential header lines if they exist
                continue

            # Check if we already have a hit recorded for this query ID
            if query_id not in best_hits:
                best_hits[query_id] = (evalue, bitscore, line)
            else:
                best_evalue, best_bitscore, _ = best_hits[query_id]
                
                # Update if the current line has a better e-value, 
                # or the same e-value but a higher bitscore
                if evalue < best_evalue or (evalue == best_evalue and bitscore > best_bitscore):
                    best_hits[query_id] = (evalue, bitscore, line)

    print(f"Writing best hits to {output_file}...")
    
    # Write the filtered lines to the output file
    with open(output_file, 'w') as outfile:
        for query_id in best_hits:
            outfile.write(best_hits[query_id][2])
            
    print(f"Done! Filtered down to {len(best_hits)} unique query genes.")

if __name__ == "__main__":
    # You can modify these file names if they differ
    input_filename = "G:\\VM-WGCNA\\WGCNA_GO-KEGG\\PURPLE_MODULE\\blastp_results.txt"
    output_filename = "G:\\VM-WGCNA\\WGCNA_GO-KEGG\\PURPLE_MODULE\\blastp_results_filtered.txt"
    
    filter_best_orthologs(input_filename, output_filename)