from Bio import SeqIO

def calculate_gc_content(sequence):
    """Calculate the GC content of a DNA sequence."""
    g_count = sequence.count("G")
    c_count = sequence.count("C")
    total_bases = len(sequence)
    if total_bases == 0:
        return 0
    return ((g_count + c_count) / total_bases) * 100

def gc_content_from_fasta(fasta_file, output_file):
    """Calculate GC content for each sequence in a FASTA file and save it in tabular format."""
    with open(output_file, "w") as out_file:
        out_file.write("GeneID\tGC_Content (%)\n")
        for record in SeqIO.parse(fasta_file, "fasta"):
            gene_id = record.id
            sequence = str(record.seq).upper()
            gc_content = calculate_gc_content(sequence)
            out_file.write(f"{gene_id}\t{gc_content:.2f}\n")

if __name__ == "__main__":
    # Specify input and output file paths
    fasta_file = "G:\Bioinformatics_Data\mango_sugar_transproter\Micro Synteny\Mdomestica\Mdomestica_transcript.fa"  # Replace with your input FASTA file path
    output_file = "G:\Bioinformatics_Data\mango_sugar_transproter\Micro Synteny\Mdomestica\Mdomestica_GC.txt"  # Replace with your desired output file path

    try:
        gc_content_from_fasta(fasta_file, output_file)
        print(f"GC content calculated and saved to {output_file}")
    except Exception as e:
        print(f"An error occurred: {e}")
