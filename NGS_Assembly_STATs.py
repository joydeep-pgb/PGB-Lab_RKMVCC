from Bio import SeqIO
import statistics

def calculate_gc_content(sequence):
    """
    Calculate the GC content percentage of a given sequence.
    
    Args:
    sequence (str): DNA sequence
    
    Returns:
    float: GC content as a percentage
    """
    g_count = sequence.count('G')
    c_count = sequence.count('C')
    gc_content = 100.0 * (g_count + c_count) / len(sequence)
    return gc_content

def analyze_fasta(file_path):
    """
    Analyze the gene length and GC content for each sequence in a FASTA file.
    
    Args:
    file_path (str): Path to the FASTA file
    
    Returns:
    dict: A dictionary where keys are sequence IDs and values are tuples of (length, GC content)
    list: A list of sequence lengths
    list: A list of GC contents
    """
    results = {}
    lengths = []
    gc_contents = []
    
    for record in SeqIO.parse(file_path, "fasta"):
        seq_id = record.id
        sequence = str(record.seq)
        seq_length = len(sequence)
        gc_content = calculate_gc_content(sequence)
        
        results[seq_id] = (seq_length, gc_content)
        lengths.append(seq_length)
        gc_contents.append(gc_content)
    
    return results, lengths, gc_contents

def calculate_n50(lengths):
    """
    Calculate the N50 of the sequence lengths, considering the total genome length.
    
    Args:
    lengths (list): List of sequence lengths
    
    Returns:
    int: The N50 value
    """
    sorted_lengths = sorted(lengths, reverse=True)
    cumulative_length = 0
    total_length = sum(sorted_lengths)
    
    for length in sorted_lengths:
        cumulative_length += length
        if cumulative_length >= total_length / 2:
            return length
    return 0  # In case the lengths list is empty

def find_second_lowest(lengths):
    """
    Find the second lowest sequence length.
    
    Args:
    lengths (list): List of sequence lengths
    
    Returns:
    int: The second lowest sequence length, or None if not applicable
    """
    unique_lengths = sorted(set(lengths))
    if len(unique_lengths) >= 2:
        return unique_lengths[1]
    return None  # Return None if there's no second unique length

def save_results_to_file(results, lengths, gc_contents, n50, second_lowest, output_file):
    """
    Save the results of the gene length and GC content analysis, including statistics, N50, and second lowest length, to a text file.
    
    Args:
    results (dict): Dictionary with sequence IDs as keys and tuples of (length, GC content) as values
    lengths (list): List of sequence lengths
    gc_contents (list): List of GC contents
    n50 (int): The N50 value of the sequence lengths
    second_lowest (int): The second lowest sequence length
    output_file (str): Path to the output text file
    """
    with open(output_file, 'w') as f:
        # Save per-sequence details
        for seq_id, (length, gc_content) in results.items():
            f.write(f"Sequence ID: {seq_id}\n")
            f.write(f"Length: {length}\n")
            f.write(f"GC Content: {gc_content:.2f}%\n\n")
        
        # Calculate and save summary statistics
        total_sequences = len(lengths)
        min_length = min(lengths)
        max_length = max(lengths)
        mean_length = statistics.mean(lengths)
        mean_gc_content = statistics.mean(gc_contents)
        
        f.write("Summary Statistics:\n")
        f.write(f"Total Number of Sequences: {total_sequences}\n")
        f.write(f"Minimum Length: {min_length}\n")
        f.write(f"Second Lowest Length: {second_lowest}\n")
        f.write(f"Maximum Length: {max_length}\n")
        f.write(f"Mean Length: {mean_length:.2f}\n")
        f.write(f"Mean GC Content: {mean_gc_content:.2f}%\n")
        f.write(f"N50 Contig Length: {n50}\n")

def main():
    # Path to your FASTA file
    fasta_file = "RF_Merge.fasta"
    # Path to the output text file
    output_file = "RF_STATs.txt"
    
    # Analyze the FASTA file
    results, lengths, gc_contents = analyze_fasta(fasta_file)
    
    # Calculate N50
    n50 = calculate_n50(lengths)
    
    # Find second lowest length
    second_lowest = find_second_lowest(lengths)
    
    # Save results and statistics to a text file
    save_results_to_file(results, lengths, gc_contents, n50, second_lowest, output_file)

if __name__ == "__main__":
    main()
