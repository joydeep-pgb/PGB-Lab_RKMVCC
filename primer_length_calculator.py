def calculate_primer_lengths(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    pair_count = 0
    total_nt_count = 0

    for i in range(0, len(lines), 2):  # Iterate through pairs of lines
        if i < len(lines) and "F:" in lines[i]:  # Check if it's a forward primer line
            forward_primer_line = lines[i].strip()
            if i + 1 < len(lines) and "R:" in lines[i + 1]:  # Check if the next line is a reverse primer
                reverse_primer_line = lines[i + 1].strip()

                # Extract sequences from the lines
                forward_primer = forward_primer_line.split(": ")[1]
                reverse_primer = reverse_primer_line.split(": ")[1]

                # Calculate lengths
                forward_length = len(forward_primer)
                reverse_length = len(reverse_primer)

                # Update the total nt count
                total_nt_count += forward_length + reverse_length

                # Output the results
                print(f"F: {forward_primer} (Length: {forward_length})")
                print(f"R: {reverse_primer} (Length: {reverse_length})")
                print("-" * 40)
                
                pair_count += 1
            else:
                print(f"Error: Expected reverse primer for forward primer on line {i + 1}.")
                break
        else:
            print(f"Error: Expected forward primer on line {i}.")
            break

    print(f"Total primer pairs processed: {pair_count}")
    print(f"Total nucleotide (nt) count: {total_nt_count}")

# Input text file with primers
filename = 'primers.txt'

# Call the function to calculate and display primer lengths
calculate_primer_lengths(filename)
