import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.stats import zscore
import numpy as np

# Read the data from the CSV file
df = pd.read_csv('Ortho_fig.txt', sep="\t")

# Display the first few rows to verify the data
print(df.head())

# Calculate the Z-score for the 'Value' column
df['Z_Value'] = zscore(df['Value'])

# Apply log10 transformation to the 'Value' column
df['Log10_Value'] = np.log10(df['Value'])

# Verify the Log10 transformation
print(df[['Category', 'Value', 'Log10_Value']].head())
print(df['Log10_Value'].describe())

# Specify the order of categories
category_order = ["Prunus dulcis", "Vitis vinifera", "Malus domestica", "Actinidia chinensis",
    "Prunus avium", "Juglans regia", "Citrus clementina", "Populus trichocarpa",
    "Chenopodium quinoa", "Vigna radiata", "Prunus persica", "Medicago truncatula",
    "Corylus avellana", "Solanum tuberosum", "Ipomoea triloba", "Manihot esculenta",
    "Cannabis sativa", "Gossypium raimondii", "Brassica rapa", "Glycine max",
    "Arabidopsis thaliana", "Triticum aestivum", "Sorghum bicolor", "Oryza sativa",
    "Zea mays"]

# Creating the plot
plt.figure(figsize=(14, 8))  # Adjust figure size for better spacing

# Jitter plot with hue assigned to 'Category' and legend suppressed
sns.stripplot(x='Category', y='Log10_Value', data=df, hue='Category', jitter=0.3, palette="Set2", size=7, alpha=0.5, legend=False, order=category_order)

# Add median line spanning the entire width of each category
for i, category in enumerate(category_order):
    if category in df['Category'].unique():
        median = df[df['Category'] == category]['Log10_Value'].median()
        plt.plot([i - 0.4, i + 0.4], [median, median], color='black', linewidth=3)

    # Count number of values in each category
    count = df[df['Category'] == category].shape[0]
    plt.text(i, df['Log10_Value'].max() + 0.1, f'{count}', ha='center', va='bottom')

# Adjust x-axis labels rotation and alignment
plt.xticks(rotation=45, ha='right')

# Remove top and right spines (axes lines)
sns.despine()

plt.title('')
plt.xlabel('Species')
plt.ylabel('Log10(Transcript Length)')
plt.tight_layout()  # Adjust layout to prevent clipping of labels
plt.show()
