import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import MultipleLocator
from matplotlib.patches import Patch

# ---------------------------
# 0️⃣ Theming
# ---------------------------
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)

# ---------------------------
# 1️⃣ Data
# ---------------------------
data = [
    ["snRNPs", 7, 12, "Drought"],
    ["SFs", 5, 4, "Drought"],
    ["Sm & Sm-like", 0, 1, "Drought"],
    ["NTC/NTR", 1, 0, "Drought"],
    ["Auxiliary proteins", 1, 3, "Drought"],

    ["snRNPs", 18, 3, "Salt"],
    ["SFs", 5, 3, "Salt"],
    ["Sm & Sm-like", 2, 0, "Salt"],
    ["NTC/NTR", 0, 1, "Salt"],
    ["Auxiliary proteins", 4, 1, "Salt"],

    ["snRNPs", 5, 9, "Heat"],
    ["SFs", 2, 5, "Heat"],
    ["Sm & Sm-like", 0, 0, "Heat"],
    ["NTC/NTR", 0, 0, "Heat"],
    ["Auxiliary proteins", 3, 0, "Heat"]
]

df = pd.DataFrame(data, columns=["Group", "UP", "Down", "Stress"])

# ---------------------------
# 2️⃣ Order
# ---------------------------
groups = ["snRNPs", "Sm & Sm-like", "NTC/NTR", "SFs", "Auxiliary proteins"]
stresses = ["Drought", "Salt", "Heat"]

# ---------------------------
# 3️⃣ Colors
# ---------------------------
colors = {
    ("Drought", "UP"):   "#F1C40F",
    ("Drought", "Down"): "#F9E79F",
    ("Salt", "UP"):      "#52BE80",
    ("Salt", "Down"):    "#A9DFBF",
    ("Heat", "UP"):      "#5DADE2",
    ("Heat", "Down"):    "#AED6F1",
}

# ---------------------------
# 4️⃣ Plot
# ---------------------------
fig, ax = plt.subplots(figsize=(6, 8))

y = np.arange(len(groups))
bar_height = 0.25

for i, stress in enumerate(stresses):
    subset = df[df["Stress"] == stress].set_index("Group").reindex(groups)

    up = subset["UP"].values
    down = subset["Down"].values
    total = up + down

    pos = y + (i - 1) * bar_height

    # UP
    ax.barh(
        pos, up, bar_height,
        color=colors[(stress, "UP")],
        linewidth=0.5
    )

    # Down (stacked)
    ax.barh(
        pos, down, bar_height,
        left=up,
        color=colors[(stress, "Down")],
        linewidth=0.5
    )

    # ---------------------------
    # 5️⃣ Annotations
    # ---------------------------
    for y_pos, u, d, t in zip(pos, up, down, total):

        # UP label (centered)
        if u > 0:
            ax.text(
                u / 2, y_pos, str(int(u)),
                va="center", ha="center",
                fontsize=9, color="black"
            )

        # Down label (centered)
        if d > 0:
            ax.text(
                u + d / 2, y_pos, str(int(d)),
                va="center", ha="center",
                fontsize=9, color="black"
            )

        # Total label (end of bar)
        if t > 0:
            ax.text(
                t + 0.3, y_pos, str(int(t)),
                va="center", ha="left",
                fontsize=9, color="black"
            )

# ---------------------------
# 6️⃣ Axes formatting
# ---------------------------
ax.set_yticks(y)
ax.set_yticklabels(groups, fontsize=14)
ax.set_xlabel("Differentially expressed spliceosome", fontsize=14)

ax.xaxis.set_major_locator(MultipleLocator(5))
#ax.xaxis.set_minor_locator(MultipleLocator(1))

ax.tick_params(axis="x", labelsize=14)
ax.tick_params(axis="y", length=0)

ax.invert_yaxis()
ax.spines["left"].set_visible(True)

# ---------------------------
# 7️⃣ Legend
# ---------------------------
legend_elements = [
    Patch(facecolor=colors[("Drought","UP")],  edgecolor="black", label="UP"),
    Patch(facecolor=colors[("Drought","Down")],edgecolor="black", label="Down"),
    Patch(facecolor=colors[("Salt","UP")],     edgecolor="black", label="UP"),
    Patch(facecolor=colors[("Salt","Down")],   edgecolor="black", label="Down"),
    Patch(facecolor=colors[("Heat","UP")],     edgecolor="black", label="UP"),
    Patch(facecolor=colors[("Heat","Down")],   edgecolor="black", label="Down"),
]

ax.legend(
    handles=legend_elements,
    loc="lower right",
    frameon=False,
    fancybox=False,      # ⬅️ critical: removes rounded corners
    handlelength=1.0,    # ⬅️ square width
    handleheight=1.0,    # ⬅️ square height
    borderpad=0.6,
    labelspacing=0.6
)

plt.tight_layout()
# Save as PDF
#plt.savefig("Enzyme_AI_NI_04.pdf", format="pdf", bbox_inches='tight')
plt.show()
