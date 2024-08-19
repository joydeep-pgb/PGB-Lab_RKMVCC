from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

def calculate_gc_content(sequence):
    # Convert fraction to percentage
    return gc_fraction(sequence) * 100

def main(file_path, output_file):
    with open(output_file, 'w') as outfile:
        outfile.write("Gene_ID\tGC Content (%)\n")
        
        # Using SeqIO to parse the FASTA file
        for record in SeqIO.parse(file_path, "fasta"):
            gc_content = calculate_gc_content(record.seq)
            outfile.write(f"{record.id}\t{gc_content:.2f}\n")

if __name__ == "__main__":
    file_path = "MangolncRNA/mango_lncRNA_Denovo.fasta"  # Replace with the path to your FASTA file
    output_file = "gc_content_results.txt"  # Define the output file name
    main(file_path, output_file)
    print(f"GC content results have been saved to {output_file}.")