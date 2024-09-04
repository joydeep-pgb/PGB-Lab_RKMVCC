import matplotlib.pyplot as plt
import numpy as np

# Data
labels = ['Genic', 'Intergenic', 'Intronic', 'lncNAT']
values = [58, 337, 48, 191]
colors = ['#ff9999', '#66b3ff', '#99ff99', '#ffcc99']

# Calculate the percentages
total = sum(values)
percentages = [v / total * 100 for v in values]

# Number of variables
num_vars = len(labels)

# Compute the angle for each bar
angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()

# Complete the circle
angles += angles[:1]
values += values[:1]
percentages += percentages[:1]

# Create the polar plot
fig, ax = plt.subplots(figsize=(6, 6), subplot_kw=dict(polar=True))

# Plot each bar
bars = ax.bar(angles[:-1], values[:-1], width=2 * np.pi / num_vars, color=colors, edgecolor='black', linewidth=1)

# Add labels with percentage
for bar, angle, percentage in zip(bars, angles[:-1], percentages[:-1]):
    height = bar.get_height()

    # Position the text slightly above the bar
    offset = height * 0.8
    ax.text(
        angle, height + offset, f'{percentage:.1f}%', ha='center', va='bottom', fontsize=12, color='black'
    )

# Add category labels and align them automatically
ax.set_xticks(angles[:-1])
ax.set_xticklabels(labels, fontsize=12)

for label, angle in zip(ax.get_xticklabels(), angles[:-1]):
    rotation = np.degrees(angle)
    x, y = label.get_position()
    
    # Dynamically adjust label positions to reduce overlap
    if rotation < 90 or rotation > 270:
        label.set_horizontalalignment('left')
        label.set_verticalalignment('center')
        label.set_position((x + 0.1, y))  # Shift right
    else:
        label.set_horizontalalignment('right')
        label.set_verticalalignment('center')
        label.set_position((x - 0.1, y))  # Shift left

# Hide the radial labels (the circle labels)
ax.set_yticklabels([])

# Customize the axis color and width
ax.spines['polar'].set_visible(False)
ax.spines['polar'].set_color('black')  # Set the desired color
ax.spines['polar'].set_linewidth(2)  # Set the desired width

# Set grid visibility and color
ax.grid(True, color='gray', linestyle='--', linewidth=0.5)

# Show the plot
plt.show()
