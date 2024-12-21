import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.interpolate import make_interp_spline
import matplotlib.colors as mcolors

# Given data
rf = np.array([0.4, 0.42, 0.52, 0.6, 0.62, 0.4, 0.42, 0.52, 0.6, 0.62, 
               0.4, 0.42, 0.52, 0.6, 0.62, 0.4, 0.42, 0.52, 0.6, 0.62, 
               0.4, 0.42, 0.52, 0.6, 0.62, 0.4, 0.42, 0.52, 0.6, 0.62])

au = np.array([0, 1, 173.7, 1, 0, 0, 1, 250, 1, 0, 0, 1, 350, 1, 0, 
               0, 1, 250.9, 1, 0, 0, 1, 270.5, 1, 0, 0, 1, 190.4, 1, 0])

track = np.array([10, 10, 10, 10, 10, 35, 35, 35, 35, 35, 60, 60, 60, 60, 60, 
                  130, 130, 130, 130, 130, 155, 155, 155, 155, 155, 180, 180, 180, 180, 180])


# Define colors for each track with hex codes and transparency
track_colors = {
    10: (*mcolors.hex2color('#0077b6'), 0.3),    
    35: (*mcolors.hex2color('#0077b6'), 0.3),    
    60: (*mcolors.hex2color('#0077b6'), 0.3),    
    130: (*mcolors.hex2color('#e9d8a6'), 0.8),    
    155: (*mcolors.hex2color('#ee9b00'), 0.8),    
    180: (*mcolors.hex2color('#9b2226'), 0.8)   
}

# Create unique tracks in reverse order
tracks = np.unique(track)[::-1]

# Prepare the 3D plot
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

for t in tracks:
    # Filter data for each track
    rf_t = rf[track == t]
    au_t = au[track == t]

    # Interpolate for smooth curve
    x_new = np.linspace(rf_t.min(), rf_t.max(), 300)
    spline = make_interp_spline(rf_t, au_t, k=3)  # Cubic spline interpolation
    z_smooth = spline(x_new)

    # Ensure no negative z-axis values
    z_smooth = np.maximum(z_smooth, 0)

    # Create array for track to match dimensions
    track_array = np.full_like(x_new, t)

    # Plot the 3D line
    ax.plot(xs=track_array, ys=x_new, zs=z_smooth, label=f"Track {t}", linewidth=2)

    # Fill under the curve using Poly3DCollection for 3D fill
    verts = [(t, x, z) for x, z in zip(x_new, z_smooth)] + [(t, x_new[-1], 0), (t, x_new[0], 0)]
    poly = Poly3DCollection([verts], color=track_colors[t], alpha=track_colors[t][-1])
    ax.add_collection3d(poly)

# Label axes
ax.set_xlabel("Track Distance (mm)", fontsize=12, labelpad=10)
ax.set_ylabel("Rf", fontsize=12, labelpad=10)
ax.set_zlabel("AU", fontsize=12, labelpad=10)

# Configure background settings
ax.grid(False)
ax.xaxis.pane.set_edgecolor('black')
ax.yaxis.pane.set_edgecolor('black')
ax.xaxis.pane.set_linewidth(1.5)
ax.yaxis.pane.set_linewidth(1.5)
ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False

# Configure axis colors and thickness
for axis in ['x', 'y', 'z']:
    ax_line = getattr(ax, f'{axis}axis').line
    ax_line.set_color('#274c77')
    ax_line.set_linewidth(2)

# Display the plot
plt.show()
