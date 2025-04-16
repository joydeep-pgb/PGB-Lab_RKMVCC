import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statannotations.Annotator import Annotator
from scipy.stats import ttest_ind

# Create the DataFrame for all genes
data_frames = []

# MiSFP7
data_frames.append(pd.DataFrame({
    'Group': ['LF', 'LF', 'LF', 'GF', 'GF', 'GF', 'RF', 'RF', 'RF'],
    'FC': [2.56, 1.88, 2.69, 4.41, 3.94, 3.78, 6.85, 5.98, 7.42],
    'Gene': 'MiSFP7'
}))

# MivGT2
data_frames.append(pd.DataFrame({
    'Group': ['LF', 'LF', 'LF', 'GF', 'GF', 'GF', 'RF', 'RF', 'RF'],
    'FC': [0.15, 1.75, 0.85, 1.52, 2.20, 1.25, 4.66, 5.34, 6.55],
    'Gene': 'MivGT2'
}))

# MivGT3
data_frames.append(pd.DataFrame({
    'Group': ['LF', 'LF', 'LF', 'GF', 'GF', 'GF', 'RF', 'RF', 'RF'],
    'FC': [2.54, 1.78, 1.52, 4.36, 4.58, 5.25, 4.84, 4.42, 5.94],
    'Gene': 'MivGT3'
}))

# MiSTP10
data_frames.append(pd.DataFrame({
    'Group': ['LF', 'LF', 'LF', 'GF', 'GF', 'GF', 'RF', 'RF', 'RF'],
    'FC': [2.36, 3.94, 2.66, 1.28, 0.95, 1.66, 3.25, 2.48, 3.55],
    'Gene': 'MiSTP10'
}))

# MiSWEET1
data_frames.append(pd.DataFrame({
    'Group': ['LF', 'LF', 'LF', 'GF', 'GF', 'GF', 'RF', 'RF', 'RF'],
    'FC': [2.04, 1.58, 1.32, 6.88, 4.82, 6.36, 1.82, 2.64, 3.26],
    'Gene': 'MiSWEET1'
}))

# MiSWEET14
data_frames.append(pd.DataFrame({
    'Group': ['LF', 'LF', 'LF', 'GF', 'GF', 'GF', 'RF', 'RF', 'RF'],
    'FC': [3.18, 4.43, 4.28, 3.86, 3.94, 6.5, 3.65, 2.36, 2.82],
    'Gene': 'MiSWEET14'
}))

# MiPMT2
data_frames.append(pd.DataFrame({
    'Group': ['LF', 'LF', 'LF', 'GF', 'GF', 'GF', 'RF', 'RF', 'RF'],
    'FC': [-0.56, -1.35, 0.66, 1.64, 2.35, 3.27, 3.56, 4.43, 2.78],
    'Gene': 'MiPMT2'
}))

# MiPMT3
data_frames.append(pd.DataFrame({
    'Group': ['LF', 'LF', 'LF', 'GF', 'GF', 'GF', 'RF', 'RF', 'RF'],
    'FC': [-1.16, -1.44, -2.68, 0.58, 0.86, 1.83, 2.84, 3.45, 3.22],
    'Gene': 'MiPMT3'
}))

# MiSUC4
data_frames.append(pd.DataFrame({
    'Group': ['LF', 'LF', 'LF', 'GF', 'GF', 'GF', 'RF', 'RF', 'RF'],
    'FC': [-0.58, 0.22, 0.62, 2.64, 3.42, 2.88, 4.21, 3.79, 3.85],
    'Gene': 'MiSUC4'
}))

# MiINT4
data_frames.append(pd.DataFrame({
    'Group': ['LF', 'LF', 'LF', 'GF', 'GF', 'GF', 'RF', 'RF', 'RF'],
    'FC': [1.64, 0.88, 1.32, 4.46, 2.76, 2.94, 2.84, 5.64, 5.08],
    'Gene': 'MiINT4'
}))

# MiSUC1
data_frames.append(pd.DataFrame({
    'Group': ['LF', 'LF', 'LF', 'GF', 'GF', 'GF', 'RF', 'RF', 'RF'],
    'FC': [1.25, 0.84, 2.64, 4.62, 3.55, 5.41, 5.89, 3.86, 3.17],
    'Gene': 'MiSUC1'
}))

# MiSUC3
data_frames.append(pd.DataFrame({
    'Group': ['LF', 'LF', 'LF', 'GF', 'GF', 'GF', 'RF', 'RF', 'RF'],
    'FC': [2.37, 1.8, 1.18, 4.04, 4.48, 2.94, 5.37, 3.26, 5.05],
    'Gene': 'MiSUC3'
}))

# MiSUC5
data_frames.append(pd.DataFrame({
    'Group': ['LF', 'LF', 'LF', 'GF', 'GF', 'GF', 'RF', 'RF', 'RF'],
    'FC': [2.49, 3.63, 1.95, 3.88, 4.47, 5.32, 1.43, 1.18, 1.02],
    'Gene': 'MiSUC5'
}))

# MiSWEET10
data_frames.append(pd.DataFrame({
    'Group': ['LF', 'LF', 'LF', 'GF', 'GF', 'GF', 'RF', 'RF', 'RF'],
    'FC': [1.61, 4.83, 2.49, 6.65, 5.85, 7.36, 1.15, 2.06, 0.48],
    'Gene': 'MiSWEET10'
}))

# Combine all DataFrames
df_all = pd.concat(data_frames, ignore_index=True)

# Set categorical order for Groups and Genes
df_all['Group'] = pd.Categorical(df_all['Group'], categories=['LF', 'GF', 'RF'], ordered=True)
genes_order = ["MiSFP7", "MivGT2", "MivGT3", "MiSTP10", "MiSWEET1", "MiSWEET14",
               "MiPMT2", "MiPMT3", "MiSUC4", "MiINT4", "MiSUC1", "MiSUC3",
               "MiSUC5", "MiSWEET10"]
df_all['Gene'] = pd.Categorical(df_all['Gene'], categories=genes_order, ordered=True)
df_all = df_all.sort_values('Gene')

# Setup figure and axes
n_rows = 2
n_cols = 7
fig, axes = plt.subplots(n_rows, n_cols, figsize=(13, 15))
axes = axes.flatten()

# Define groups and colors
groups = ['LF', 'GF', 'RF']
palette = {'LF': '#5aa39e', 'GF': '#fddd5d', 'RF': '#f4f3f1'}

# Create each subplot
for i, (gene, ax) in enumerate(zip(df_all['Gene'].cat.categories, axes)):
    data_sub = df_all[df_all['Gene'] == gene]

    # Boxplot
    sns.boxplot(x='Group', y='FC', data=data_sub, ax=ax, palette=palette, order=groups,
                width=0.6, linewidth=1, fliersize=0,
                boxprops=dict(edgecolor='black', linewidth=1),
                whiskerprops=dict(color='black', linewidth=1),
                capprops=dict(color='black', linewidth=1),
                medianprops=dict(color='black', linewidth=1))

    # Jitter points
    sns.stripplot(x='Group', y='FC', data=data_sub, ax=ax, order=groups,
                  jitter=0.2, color='black', edgecolor='white', linewidth=1.3,
                  facecolor='black', size=8, alpha=1)

    # Statistical annotations
    pairs = [('LF', 'GF'), ('GF', 'RF'), ('LF', 'RF')]
    annotator = Annotator(ax, pairs, data=data_sub, x='Group', y='FC', order=groups)
    annotator.configure(test='t-test_ind', text_format='star',
                       comparisons_correction=None, loc='outside',
                       fontsize=12, line_height=0.03, line_width=1)
    annotator.apply_and_annotate()

    # Subplot aesthetics
    ax.set_xlabel(gene, fontsize=12, labelpad=10)
    ax.set_title('')
    ax.set_ylabel('Relative Fold Change' if i % n_cols == 0 else '', fontsize=14)
    ax.tick_params(axis='x', labelsize=12, color='#4a4e69', length=6, width=1.5)
    ax.tick_params(axis='y', labelsize=12, color='#4a4e69', length=6, width=1.5)

    # Spine adjustments
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    for spine in ['bottom', 'left']:
        ax.spines[spine].set_color('#4a4e69')
        ax.spines[spine].set_linewidth(1.5)

# Hide empty subplots if any
for j in range(len(df_all['Gene'].cat.categories), len(axes)):
    axes[j].axis('off')

plt.tight_layout()

# Save as PDF
plt.savefig("Gene_FC_Plots.pdf", format="pdf", bbox_inches='tight')

# Display the plot
plt.show()
