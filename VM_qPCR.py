import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from statannotations.Annotator import Annotator
from matplotlib.ticker import MultipleLocator

# --- Set Global Font to Arial ---
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

# Create DataFrame
data = {
    "lncRNA_FC": [3.64, 5.82, 5.28, 1.16, 0.84, 1.48],
    "mRNA_FC":   [5.66, 4.80, 6.46, 3.22, 1.08, 2.44],
    "Group":     ["VM-R", "VM-R", "VM-R", "VM-S", "VM-S", "VM-S"]
}
df = pd.DataFrame(data)

# Melt to long format
df_long = pd.melt(df, id_vars=["Group"], value_vars=["lncRNA_FC", "mRNA_FC"],
                  var_name="GeneType", value_name="FC")
df_long["GeneType"] = df_long["GeneType"].replace({
    "lncRNA_FC": "VmLnc.41055.1",
    "mRNA_FC":   "ANK"
})

# Plot styling
sns.set(style="white", font_scale=1.2)
palette = {"VM-R": "#b9d14c", "VM-S": "#358ea5"}
order = ["VmLnc.41055.1", "ANK"]
group_order = ["VM-R", "VM-S"]

fig, ax = plt.subplots(figsize=(3, 5.5))

# Boxplot
sns.boxplot(x='GeneType', y='FC', data=df_long, ax=ax, hue='Group',
            palette=palette, order=order, width=0.65,
            linewidth=1.2, fliersize=0,
            boxprops=dict(edgecolor='black', linewidth=1.2),
            whiskerprops=dict(color='black', linewidth=1.2),
            capprops=dict(color='black', linewidth=1.2),
            medianprops=dict(color='black', linewidth=1.2))

# Remove the legend completely
if ax.get_legend():
    ax.get_legend().remove()

# White circle mean markers
means = df_long.groupby(['GeneType', 'Group'], observed=False)['FC'].mean().reset_index()
for i, gene in enumerate(order):
    for j, group in enumerate(group_order):
        val = means[(means['GeneType'] == gene) & (means['Group'] == group)]['FC'].values
        if val.size > 0:
            xpos = i - 0.165 if group == group_order[0] else i + 0.165
            ax.plot(xpos, val[0], marker='o', markersize=8,
                    markeredgecolor='black', markerfacecolor='white',
                    markeredgewidth=1.3, zorder=10)

# Statistical annotation
# NOTE: If you decide to add the other pair back to test both, 
# the new configure settings will automatically hide it if it returns "ns".
pairs = [(("VmLnc.41055.1", "VM-R"), ("VmLnc.41055.1", "VM-S")), 
          (("ANK", "VM-R"), ("ANK", "VM-S"))]

annotator = Annotator(ax, pairs, data=df_long, x='GeneType', y='FC',
                      hue='Group', order=order, hue_order=group_order)

# Added hide_non_significant=True to automatically skip ns brackets
annotator.configure(test='t-test_ind', text_format='star',
                    comparisons_correction=None, loc='outside',
                    fontsize=14, line_height=0.03, line_width=1.2,
                    hide_non_significant=True)

annotator.apply_and_annotate()

# X-axis: 4 ticks at each box position
tick_positions = [0 - 0.165, 0 + 0.165, 1 - 0.165, 1 + 0.165]
tick_labels    = ["VM-R", "VM-S", "VM-R", "VM-S"]
ax.set_xticks(tick_positions)
ax.set_xticklabels(tick_labels, fontsize=10.5)

# Gene names centred under each pair
for i, gene in enumerate(order):
    ax.text(i, -0.13, gene, ha='center', va='top', fontsize=12,
            fontweight='bold', transform=ax.get_xaxis_transform())

# Axis labels
ax.set_xlabel("")
ax.set_ylabel("", fontsize=13)
ax.set_title("")
ax.set_ylim(0, None)

ax.set_ylim(0, None)
ax.yaxis.set_major_locator(MultipleLocator(2))

# Spine customization
for spine in ['top', 'right']:
    ax.spines[spine].set_visible(False)
for spine in ['left', 'bottom']:
    ax.spines[spine].set_color('#4a4e69')
    ax.spines[spine].set_linewidth(1.2)

# Force inward tick marks on both axes
ax.tick_params(axis='x', which='both', bottom=True,
               direction='out', length=3.5, width=1.2, color='#4a4e69', labelsize=10.5)
ax.tick_params(axis='y', which='both', left=True,
               direction='out', length=3.5, width=1.2, color='#4a4e69', labelsize=11)

# Editable text for Illustrator
plt.rcParams['pdf.fonttype'] = 42
plt.tight_layout()
plt.subplots_adjust(bottom=0.15)
plt.show()