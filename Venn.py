import matplotlib.pyplot as plt
from matplotlib_venn import venn3

# Function to read gene IDs from a file and return a set
def read_genes(file_path):
    with open(file_path, 'r') as file:
        genes = set(line.strip() for line in file)
    return genes

# Load gene sets from the files
a = read_genes('LFvsGF.txt')
b = read_genes('GFvsRF.txt')
c = read_genes('LFvsRF.txt')

# Create the Venn diagram
plt.figure(figsize=(8, 8))
venn = venn3([a, b, c], ('a', 'b', 'c'))

# Set the title for the diagram
plt.title("Proportional Venn Diagram of Gene Sets")

# Annotating the intersections with labels and counts
if venn.get_label_by_id('100'):
    venn.get_label_by_id('100').set_text(f'a\n({len(a - b - c)})')
if venn.get_label_by_id('010'):
    venn.get_label_by_id('010').set_text(f'b\n({len(b - a - c)})')
if venn.get_label_by_id('001'):
    venn.get_label_by_id('001').set_text(f'c\n({len(c - a - b)})')
if venn.get_label_by_id('110'):
    venn.get_label_by_id('110').set_text(f'ab\n({len(a & b - c)})')
if venn.get_label_by_id('101'):
    venn.get_label_by_id('101').set_text(f'ac\n({len(a & c - b)})')
if venn.get_label_by_id('011'):
    venn.get_label_by_id('011').set_text(f'bc\n({len(b & c - a)})')
if venn.get_label_by_id('111'):
    venn.get_label_by_id('111').set_text(f'abc\n({len(a & b & c)})')

# Display the diagram
plt.show()
