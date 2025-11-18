def calculate_primer_lengths(filename):
    with open(filename, 'r') as file:
        lines = [line.strip() for line in file.readlines() if line.strip()]

    triple_count = 0
    total_nt_count = 0
    i = 0

    while i < len(lines):
        # Expect F line
        if not lines[i].startswith("F:"):
            print(f"Error: Expected F: on line {i+1}")
            break
        
        forward = lines[i].split(": ")[1]
        # Expect R1 line
        if i+1 >= len(lines) or not lines[i+1].startswith("R1:"):
            print(f"Error: Expected R1: after F: on line {i+2}")
            break
        
        reverse1 = lines[i+1].split(": ")[1]
        # Expect R2 line
        if i+2 >= len(lines) or not lines[i+2].startswith("R2:"):
            print(f"Error: Expected R2: after R1: on line {i+3}")
            break
        
        reverse2 = lines[i+2].split(": ")[1]

        # Calculate lengths
        f_len = len(forward)
        r1_len = len(reverse1)
        r2_len = len(reverse2)

        # Update totals
        total_nt_count += (f_len + r1_len + r2_len)
        triple_count += 1

        # Print results
        print(f"F:  {forward} (Length: {f_len})")
        print(f"R1: {reverse1} (Length: {r1_len})")
        print(f"R2: {reverse2} (Length: {r2_len})")
        print("-" * 50)

        # Move to next block of 3 lines
        i += 3

    print(f"Total primer sets processed: {triple_count}")
    print(f"Total nucleotide count: {total_nt_count}")


# Run
filename = 'AS_Primer.txt'
calculate_primer_lengths(filename)
