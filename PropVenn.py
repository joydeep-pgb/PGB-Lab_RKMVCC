import matplotlib.pyplot as plt
from matplotlib_venn import venn2

# Your data
A = 1682
B = 2512
A_inter_B = 221

# Calculate individual sections
A_only = A - A_inter_B
B_only = B - A_inter_B

# Create figure with better proportions
plt.figure(figsize=(10, 8))

# Create the Venn diagram (subsets order: (A_only, B_only, A∩B))
venn = venn2(subsets=(A_only, B_only, A_inter_B),
             set_labels=('Set A', 'Set B'),
             set_colors=('#4e79a7', '#f28e2b'),
             alpha=0.7)

# Customize labels with formatting
venn.get_label_by_id('10').set_text(f'{A_only:,}\n({A_only/A*100:.1f}% of A)')
venn.get_label_by_id('01').set_text(f'{B_only:,}\n({B_only/B*100:.1f}% of B)')
venn.get_label_by_id('11').set_text(f'{A_inter_B:,}\n({A_inter_B/A*100:.1f}% of A\n{A_inter_B/B*100:.1f}% of B)')

# Style improvements
plt.title(f'Proportional Venn Diagram\nTotal A: {A:,} | Total B: {B:,} | Overlap: {A_inter_B:,}', 
          fontsize=14, fontweight='bold')

# Add summary text
total_union = A_only + B_only + A_inter_B
plt.text(-0.8, -0.8, f'Total unique items (A∪B): {total_union:,}', 
         fontsize=10, bbox=dict(facecolor='white', alpha=0.8))

# Adjust layout
plt.tight_layout()
plt.show()

# Print detailed statistics
print("="*60)
print("VENN DIAGRAM STATISTICS")
print("="*60)
print(f"Set A total:           {A:,}")
print(f"Set B total:           {B:,}")
print(f"Intersection (A∩B):    {A_inter_B:,}")
print(f"A only:               {A_only:,} ({A_only/A*100:.1f}% of A)")
print(f"B only:               {B_only:,} ({B_only/B*100:.1f}% of B)")
print(f"Union (A∪B):          {total_union:,}")
print(f"Jaccard similarity:   {A_inter_B/total_union:.3f}")
print(f"Overlap coefficient:  {A_inter_B/min(A, B):.3f}")
print("="*60)