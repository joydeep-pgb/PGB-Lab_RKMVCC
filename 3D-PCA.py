import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as mcolors

# Load data from CSV file
data = pd.read_csv('PCA_Plot.csv')

# Extract PCA results and groups
X_pca = data[['PC1', 'PC2', 'PC3']].values
groups = data['group']

# Define colors and markers for each group
unique_groups = groups.unique()
colors = ['#e26313', '#2a9d8f', '#ffb703', 'green', 'purple', 'orange']  # Add more as needed
markers = ['o']  # Add more as needed

# Create dictionaries for color and marker mapping
color_dict = {group: colors[i % len(colors)] for i, group in enumerate(unique_groups)}
marker_dict = {group: markers[i % len(markers)] for i, group in enumerate(unique_groups)}

# Configuration settings
marker_size = 200  # Marker size
xlabel = 'PC1, 14.89% variation'
ylabel = 'PC2, 6.74% variation'
zlabel = 'PC3, 5.65% variation'
xlabel_size, ylabel_size, zlabel_size = 15, 15, 15
tick_size = 12
label_pad = 20
axis_color = '#274c77'
axis_thickness = 2
grid_color = '#3a7ca5'
grid_thickness = 0.5
grid_alpha = 0.5
pane_color = (*mcolors.hex2color('#dbe4ee'), 0.1)  # Background pane color with transparency

# Create a 3D plot
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Plot each group
for group in unique_groups:
    indices = groups == group
    ax.scatter(X_pca[indices, 0], X_pca[indices, 1], X_pca[indices, 2],
               c=color_dict[group], marker=marker_dict[group], label=group, s=marker_size,
               edgecolors='black', linewidths=1.5)  # Added edge color and linewidths

# Label the axes
ax.set_xlabel(xlabel, fontsize=xlabel_size, labelpad=label_pad)
ax.set_ylabel(ylabel, fontsize=ylabel_size, labelpad=label_pad)
ax.set_zlabel(zlabel, fontsize=zlabel_size, labelpad=label_pad)

# Configure tick parameters
ax.tick_params(axis='both', which='major', labelsize=tick_size)

# Configure axis colors and thickness
for axis in ['x', 'y', 'z']:
    ax_line = getattr(ax, f'{axis}axis').line
    ax_line.set_color(axis_color)
    ax_line.set_linewidth(axis_thickness)

# Set the grid color, linewidth, and alpha for all axes
ax.xaxis._axinfo['grid'].update(color=grid_color, linewidth=grid_thickness, alpha=grid_alpha)
ax.yaxis._axinfo['grid'].update(color=grid_color, linewidth=grid_thickness, alpha=grid_alpha)
ax.zaxis._axinfo['grid'].update(color=grid_color, linewidth=grid_thickness, alpha=grid_alpha)

# Make the panes transparent and set the color
ax.xaxis.set_pane_color(pane_color)
ax.yaxis.set_pane_color(pane_color)
ax.zaxis.set_pane_color(pane_color)

# Add a legend
ax.legend(title="Group", loc='upper right', fontsize=12)

# Show the plot
plt.show()
