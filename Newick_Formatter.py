from ete3 import Tree

# Step 1: Read Newick tree from input file
with open(r"G:\Bioinformatics_Data\mango_sugar_transproter\mango_TA4_transcriptome\Phylogenetic tree\Final DATA\NEW\ST_New.txt", "r") as infile:
    newick_str = infile.read().strip()

# Step 2: Load the tree
tree = Tree(newick_str, format=1)

# Step 3: Add '.....' before and after each leaf (tip) name
for leaf in tree.iter_leaves():
    leaf.name = f".....{leaf.name}....."

# Step 4: Write the modified tree to a new file
with open(r"G:\Bioinformatics_Data\mango_sugar_transproter\mango_TA4_transcriptome\Phylogenetic tree\Final DATA\NEW\ST_New_Format.txt", "w") as outfile:
    outfile.write(tree.write(format=1))

print("Modified tree saved to 'modified_tree.txt'")
