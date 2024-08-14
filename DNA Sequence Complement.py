def complement_base(base):
    """Returns the complementary base."""
    if base == 'A':
        return 'T'
    elif base == 'T':
        return 'A'
    elif base == 'C':
        return 'G'
    elif base == 'G':
        return 'C'
    else:
        return base  # If the base is not A, T, C, or G, return the same base

def reverse_complement(seq):
    """Returns the reverse complement of the sequence."""
    complement_seq = [complement_base(base) for base in seq]
    reverse_complement_seq = ''.join(complement_seq[::-1])
    return reverse_complement_seq

# Get DNA sequence input from the user
dna_sequence = input("Enter the DNA sequence: ").upper()

# Calculate complement and reverse complement
complement_sequence = ''.join([complement_base(base) for base in dna_sequence])
reverse_complement_sequence = reverse_complement(dna_sequence)

# Display the results
print(f"Complement: {complement_sequence}")
print(f"Reverse complement: {reverse_complement_sequence}")
