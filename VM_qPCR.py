import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats
from statannotations.Annotator import Annotator
from matplotlib.ticker import MultipleLocator
import warnings
warnings.filterwarnings("ignore")

# --- Global font settings ---
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['pdf.fonttype'] = 42

# ─── All 11 gene pairs ───────────────────────────────────────────────────────
# Order: upper panel (6) then lower panel (5)
gene_pairs = [
    # ── Upper panel ──
    {"lncRNA": "VmLnc.42403.2", "mRNA": "bHLH",
     "lncRNA_R": [3.36, 4.70, 5.64], "lncRNA_S": [1.90, 3.60, 2.60],
     "mRNA_R":   [3.10, 5.50, 6.74], "mRNA_S":   [1.60, 2.60, 2.08]},

    {"lncRNA": "VmLnc.41325.1", "mRNA": "PR5",
     "lncRNA_R": [5.80, 8.64, 4.20], "lncRNA_S": [2.20, 3.20, 1.50],
     "mRNA_R":   [7.28, 9.10, 6.50], "mRNA_S":   [2.00, 1.20, 0.60]},

    {"lncRNA": "VmLnc.42248.9", "mRNA": "NAC",
     "lncRNA_R": [4.20, 1.90, 3.40], "lncRNA_S": [1.10, 1.80, 2.10],
     "mRNA_R":   [5.40, 4.20, 3.70], "mRNA_S":   [1.10, 0.60, 1.96]},

    {"lncRNA": "VmLnc.42513.6", "mRNA": "WRKY",
     "lncRNA_R": [3.24, 3.08, 5.94], "lncRNA_S": [3.20, 1.80, 1.36],
     "mRNA_R":   [1.80, 2.80, 3.70], "mRNA_S":   [0.30, 1.80, 0.66]},

    {"lncRNA": "VmLnc.41055.1", "mRNA": "ANK",
     "lncRNA_R": [3.64, 5.82, 5.28], "lncRNA_S": [1.16, 0.84, 1.48],
     "mRNA_R":   [5.66, 4.80, 6.46], "mRNA_S":   [3.22, 1.08, 2.44]},

    {"lncRNA": "VmLnc.28498.3", "mRNA": "TGA7",
     "lncRNA_R": [4.80, 7.60, 6.30], "lncRNA_S": [1.10, 1.70, 2.30],
     "mRNA_R":   [6.40, 7.70, 5.30], "mRNA_S":   [2.30, 1.45, 2.95]},

    # ── Lower panel ──
    {"lncRNA": "VmLnc.423.2",   "mRNA": "GDSL lipase",
     "lncRNA_R": [2.70, 1.80, 3.95], "lncRNA_S": [1.90, 2.10, 3.20],
     "mRNA_R":   [3.84, 4.62, 1.95], "mRNA_S":   [2.15, 1.20, 3.30]},

    {"lncRNA": "VmLnc.839.4",   "mRNA": "STK",
     "lncRNA_R": [1.80, 3.60, 2.30], "lncRNA_S": [0.68, 1.80, 2.50],
     "mRNA_R":   [4.40, 5.36, 2.70], "mRNA_S":   [0.90, 1.45, 1.90]},

    {"lncRNA": "VmLnc.18958.1", "mRNA": "RLK",
     "lncRNA_R": [2.30, 3.90, 5.90], "lncRNA_S": [1.60, 0.74, 0.36],
     "mRNA_R":   [6.40, 3.90, 5.72], "mRNA_S":   [1.50, 0.90, 1.70]},

    {"lncRNA": "VmLnc.22558.3", "mRNA": "CC-NB-LRR",
     "lncRNA_R": [3.90, 4.80, 2.86], "lncRNA_S": [0.60, 0.90, 1.80],
     "mRNA_R":   [2.46, 4.22, 1.90], "mRNA_S":   [0.60, 0.45, 1.85]},

    {"lncRNA": "VmLnc.33854.1", "mRNA": "WD40",
     "lncRNA_R": [1.10, 1.60, 2.30], "lncRNA_S": [0.20, 0.65, 0.70],
     "mRNA_R":   [1.84, 2.60, 2.80], "mRNA_S":   [0.30, 0.90, 0.70]},
]

upper_pairs = gene_pairs[:6]
lower_pairs = gene_pairs[6:]   # 5 pairs

palette     = {"VM-R": "#b9d14c", "VM-S": "#358ea5"}
group_order = ["VM-R", "VM-S"]


# ─── Helper: draw one subplot ────────────────────────────────────────────────
def make_subplot(ax, pair, show_ylabel=False):
    lnc  = pair["lncRNA"]
    mrna = pair["mRNA"]
    order = [lnc, mrna]

    # Build long-format DataFrame
    rows = []
    for gene, r_vals, s_vals in [(lnc,  pair["lncRNA_R"], pair["lncRNA_S"]),
                                  (mrna, pair["mRNA_R"],   pair["mRNA_S"])]:
        for v in r_vals:
            rows.append({"GeneType": gene, "Group": "VM-R", "FC": v})
        for v in s_vals:
            rows.append({"GeneType": gene, "Group": "VM-S", "FC": v})
    df = pd.DataFrame(rows)

    sns.set(style="white", font_scale=1.0)

    # Boxplot
    sns.boxplot(x='GeneType', y='FC', data=df, ax=ax, hue='Group',
                palette=palette, order=order, width=0.65, gap=0.25,
                linewidth=1.2, fliersize=0,
                boxprops=dict(edgecolor='black', linewidth=1.2),
                whiskerprops=dict(color='black', linewidth=1.2),
                capprops=dict(color='black', linewidth=1.2),
                medianprops=dict(color='black', linewidth=1.2))

    if ax.get_legend():
        ax.get_legend().remove()

    # White-circle mean markers
    means = df.groupby(['GeneType', 'Group'], observed=False)['FC'].mean().reset_index()
    for i, gene in enumerate(order):
        for group in group_order:
            val = means[(means['GeneType'] == gene) & (means['Group'] == group)]['FC'].values
            if val.size > 0:
                xpos = i - 0.165 if group == "VM-R" else i + 0.165
                ax.plot(xpos, val[0], marker='o', markersize=7,
                        markeredgecolor='black', markerfacecolor='white',
                        markeredgewidth=1.2, zorder=10)

    # Significance brackets (VM-R vs VM-S within each gene)
    pairs_stat = [((lnc,  "VM-R"), (lnc,  "VM-S")),
                  ((mrna, "VM-R"), (mrna, "VM-S"))]
    annotator = Annotator(ax, pairs_stat, data=df,
                          x='GeneType', y='FC', hue='Group',
                          order=order, hue_order=group_order)
    annotator.configure(test='t-test_ind', text_format='star',
                        comparisons_correction=None, loc='outside',
                        fontsize=11, line_height=0.025, line_width=1.0,
                        hide_non_significant=True)
    annotator.apply_and_annotate()

    # X-axis: 4 tick positions (one per box)
    tick_pos    = [0 - 0.165, 0 + 0.165, 1 - 0.165, 1 + 0.165]
    tick_labels = ["VM-R", "VM-S", "VM-R", "VM-S"]
    ax.set_xticks(tick_pos)
    ax.set_xticklabels(tick_labels, fontsize=8.5)

    # Bold gene-name labels centred under each pair
    ax.text(0, -0.11, lnc,  ha='center', va='top', fontsize=8.5,
            fontweight='bold', transform=ax.get_xaxis_transform())
    ax.text(1, -0.11, mrna, ha='center', va='top', fontsize=8.5,
            fontweight='bold', transform=ax.get_xaxis_transform())

    # Labels & limits
    ax.set_xlabel("")
    ax.set_ylabel("Relative Fold Change" if show_ylabel else "", fontsize=10)
    ax.set_title("")
    ax.set_ylim(0, None)
    ax.yaxis.set_major_locator(MultipleLocator(2))

    # Spines
    for sp in ['top', 'right']:
        ax.spines[sp].set_visible(False)
    for sp in ['left', 'bottom']:
        ax.spines[sp].set_color('#4a4e69')
        ax.spines[sp].set_linewidth(1.2)

    # Ticks
    ax.tick_params(axis='x', which='both', bottom=True,
                   direction='out', length=3.5, width=1.2,
                   color='#4a4e69', labelsize=8.5)
    ax.tick_params(axis='y', which='both', left=True,
                   direction='out', length=3.5, width=1.2,
                   color='#4a4e69', labelsize=9)


# ─── Figure layout ───────────────────────────────────────────────────────────
# 6 columns; upper row fills all 6, lower row fills 5 (left-aligned)
fig = plt.figure(figsize=(17, 12))

# Outer GridSpec: 2 rows × 6 columns
gs = gridspec.GridSpec(2, 6,
                       figure=fig,
                       hspace=0.60,   # vertical gap between rows
                       wspace=0.38)   # horizontal gap between columns

# ── Upper panel: 6 subplots ──
for i, pair in enumerate(upper_pairs):
    ax = fig.add_subplot(gs[0, i])
    make_subplot(ax, pair, show_ylabel=(i == 0))

# ── Lower panel: 5 subplots (left-aligned, same column width) ──
for i, pair in enumerate(lower_pairs):
    ax = fig.add_subplot(gs[1, i])
    make_subplot(ax, pair, show_ylabel=(i == 0))

# ─── Export ──────────────────────────────────────────────────────────────────
plt.tight_layout()
plt.show()