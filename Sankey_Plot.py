import pandas as pd
import plotly.graph_objects as go
import io

data = """LncRNA	miRNA	mRNA	GO
MSTRG.23037.3	miRNA	HD-ZIP	Response to stimulus
MSTRG.6685.20	miR166	TBL	Transferase activity
MSTRG.8053.8	miR426	DOF-ZF	DNA binding
MSTRG.8053.8	miR5658	ANK	Immune response
MSTRG.8053.8	miR5658	RIPK	Defense response
MSTRG.8053.8	miR5658	SBT5.2	Defense response
MSTRG.8053.8	miR5658	GDSL	Defense response
MSTRG.17590.1	miR866	F-box	Response to wounding
MSTRG.17113.2	miR4397	EamA	Transmembrane transport
MSTRG.34091.1	miR5035	PDI	Protein folding
MSTRG.2867.15	miR5675	TNL	Defense response"""

# Read and clean data
df = pd.read_csv(io.StringIO(data), sep='\t')
df.columns = df.columns.str.strip()
for col in df.columns:
    df[col] = df[col].astype(str).str.strip()

# Columns defining the flow stages
cols = ['LncRNA', 'miRNA', 'mRNA', 'GO']

## 1. Data Preparation for Plotly Sankey

# Initialize lists for labels and flows
all_labels = []
link_data = {'source': [], 'target': [], 'value': [], 'color': []}

# Get all unique nodes and map them to a continuous index
for col in cols:
    all_labels.extend(df[col].unique().tolist())
unique_labels = list(dict.fromkeys(all_labels))
label_to_index = {label: i for i, label in enumerate(unique_labels)}

# Prepare link data
for i in range(len(cols) - 1):
    src_col = cols[i]
    tgt_col = cols[i+1]
    
    # Group by source and target to get the flow volume (count)
    link_counts = df.groupby([src_col, tgt_col]).size().reset_index(name='count')
            
    for _, row in link_counts.iterrows():
        source_label = row[src_col]
        target_label = row[tgt_col]
        
        link_data['source'].append(label_to_index[source_label])
        link_data['target'].append(label_to_index[target_label])
        link_data['value'].append(row['count'])

## 2. Coloring Nodes and Links

# Generate a color list (Plotly uses its own default colors, but we can override)
# Using 'Pastel' palette for a soft look
import plotly.express as px
node_colors = px.colors.qualitative.Pastel

# Map each unique label to a specific color
label_to_color = {}
for i, label in enumerate(unique_labels):
    label_to_color[label] = node_colors[i % len(node_colors)]

# Create a list of node colors corresponding to the unique_labels order
node_color_list = [label_to_color[label] for label in unique_labels]

# Assign link colors based on the source node color
# Plotly link colors need to be strings (e.g., 'rgba(255, 0, 0, 0.5)')
for source_index in link_data['source']:
    source_label = unique_labels[source_index]
    hex_color = label_to_color[source_label]
    # Set link color to source node color with some transparency (alpha=0.6)
    link_data['color'].append(hex_color.replace('rgb', 'rgba').replace(')', ', 0.6)'))


## 3. Create Plotly Sankey Diagram

# Plotly automatically handles node positioning and link routing.
fig = go.Figure(data=[go.Sankey(
    # Node Configuration
    node=dict(
        pad=15,
        # Reduced thickness (Matplotlib's 0.2 width is comparable to a small thickness value in Plotly)
        thickness=10, 
        line=dict(color="black", width=0.5),
        label=unique_labels,
        color=node_color_list
    ),
    # Link Configuration
    link=dict(
        source=link_data['source'],
        target=link_data['target'],
        value=link_data['value'],
        color=link_data['color'] # Colored links!
    )
)])

fig.update_layout(
    title_text="", 
    font_size=20,
    hovermode='x' # Improves tooltip experience
)

fig.show() # This command displays the interactive plot in a compatible environment