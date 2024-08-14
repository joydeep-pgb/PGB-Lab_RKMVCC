# Script to concatenate multiple FASTA files into a single output file
import glob

def read_fasta(file):
    items = []
    index = 0
    for line in file:
        if line.startswith(">"):
            if index >= 1:
                items.append(aninstance)
            index += 1
            name = line[:-1]
            seq = ''
            aninstance = (name, seq)
        else:
            seq += line[:-1]
            aninstance = (name, seq)
            items.append(aninstance)
    return items

# Obtain directory containing single FASTA files
filepattern = input('.fasta')

# Obtain output filename
outfile = input('merge.fasta')

# Create new output file
with open(outfile, 'w') as output:
    for file in glob.glob(filepattern):
        contents = read_fasta(open(file).readlines())
        for item in contents:
            output.write(f"{item[0]}\n{item[1]}\n\n")

print("Done! Merged sequences saved in", outfile)
