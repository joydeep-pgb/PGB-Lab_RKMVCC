import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from statannotations.Annotator import Annotator

# Create DataFrame
data = {
    "lncRNA": ["VmLnc.42513.6"] * 6,
    "lncRNA_FC": [3.24, 3.08, 5.94, 3.2, 1.8, 1.36],
    "mRNA": ["WRKY"] * 6,
    "mRNA_FC": [1.8, 2.8, 3.7, 0.8, 1.5, 0.5],
    "Group": ["VM-R", "VM-R", "VM-R", "VM-S", "VM-S", "VM-S"]
}
df = pd.DataFrame(data)

# Melt to long format
df_long = pd.melt(df, id_vars=["Group"], value_vars=["lncRNA_FC", "mRNA_FC"],
                  var_name="GeneType", value_name="FC")

# Rename for clarity
df_long["GeneType"] = df_long["GeneType"].replace({
    "lncRNA_FC": "VmLnc.42513.6",
    "mRNA_FC": "WRKY"
})

# Plot styling
sns.set(style="white", font_scale=1.2)
palette = {"VM-R": "#66c2a5", "VM-S": "#fc8d62"}
order = ["VmLnc.42513.6", "WRKY"]
group_order = ["VM-R", "VM-S"]

plt.figure(figsize=(7, 5))
ax = plt.gca()

# Boxplot
sns.boxplot(x='GeneType', y='FC', data=df_long, ax=ax, hue='Group',
            palette=palette, order=order, width=0.5,
            linewidth=1, fliersize=0,
            boxprops=dict(edgecolor='black', linewidth=1),
            whiskerprops=dict(color='black', linewidth=1),
            capprops=dict(color='black', linewidth=1),
            medianprops=dict(color='black', linewidth=1))

# # Jitter points
# sns.stripplot(x='GeneType', y='FC', data=df_long, ax=ax, hue='Group',
#               order=order, dodge=True, jitter=0.2,
#               size=6, alpha=0.8, linewidth=1.2, edgecolor='white', palette='dark:.3')

# Remove duplicate legend
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[:2], labels[:2], title='Group', bbox_to_anchor=(1.05, 1), loc='upper left')

# Add mean markers manually
means = df_long.groupby(['GeneType', 'Group'], observed=False)['FC'].mean().reset_index()
for i, gene in enumerate(order):
    for j, group in enumerate(group_order):
        val = means[(means['GeneType'] == gene) & (means['Group'] == group)]['FC'].values
        if val.size > 0:
            xpos = i - 0.2 if group == group_order[0] else i + 0.2
            ax.plot(xpos, val[0], marker='o', markersize=7,
                    markeredgecolor='black', markerfacecolor='white',
                    markeredgewidth=1.3, zorder=10)

# Statistical annotation
pairs = [
    (("VmLnc.42513.6", "VM-R"), ("VmLnc.42513.6", "VM-S")),
    (("WRKY", "VM-R"), ("WRKY", "VM-S"))
]
annotator = Annotator(ax, pairs, data=df_long, x='GeneType', y='FC',
                      hue='Group', order=order, hue_order=group_order)
annotator.configure(test='t-test_ind', text_format='star',
                    comparisons_correction=None, loc='outside',
                    fontsize=13, line_height=0.03, line_width=1)
annotator.apply_and_annotate()

# Axis and title formatting
ax.set_xlabel("")
ax.set_ylabel("Fold Change (log2)", fontsize=14)
ax.set_title("Expression of VmLnc.42513.6 and WRKY", fontsize=15, weight='bold')
ax.tick_params(axis='x', labelsize=12, color='#5d637f', length=6, width=1.5)
ax.tick_params(axis='y', labelsize=12, color='#5d637f', length=6, width=1.5)

# Spine customization
for spine in ['top', 'right']:
    ax.spines[spine].set_visible(False)
for spine in ['left', 'bottom']:
    ax.spines[spine].set_color('#4a4e69')
    ax.spines[spine].set_linewidth(1.5)

# Editable text for Illustrator
plt.rcParams['pdf.fonttype'] = 42
plt.tight_layout()
plt.show()
