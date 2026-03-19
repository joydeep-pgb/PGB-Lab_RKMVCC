import pandas as pd
import requests

def convert_gprofiler_ids_to_symbols(input_csv, output_csv):
    print("Loading file...")
    df = pd.read_csv(input_csv)
    
    # 1. Extract all unique AGI codes (ATxGxxxxx) from the 'intersections' column
    all_genes = set()
    for row in df['intersections'].dropna():
        all_genes.update(row.split(','))
    
    print(f"Found {len(all_genes)} unique genes. Querying MyGene.info...")
    
    # 2. Query MyGene.info API for Arabidopsis (Taxon ID: 3702)
    url = 'https://mygene.info/v3/query'
    data = {
        'q': ','.join(all_genes),
        'scopes': 'locus_tag,tair,ensemblgene',
        'fields': 'symbol',
        'species': '3702' 
    }
    
    response = requests.post(url, data=data)
    response.raise_for_status()
    
    # 3. Create a mapping dictionary: { "AT3G15990": "AHA2", ... }
    mapping = {}
    for match in response.json():
        query_id = match.get('query')
        # If it finds a common symbol, use it. Otherwise, keep the original AT ID.
        symbol = match.get('symbol', query_id) 
        mapping[query_id] = symbol
        
    # 4. Function to replace the comma-separated IDs in the dataframe
    def replace_with_symbols(gene_string):
        if pd.isna(gene_string):
            return gene_string
        genes = gene_string.split(',')
        return ','.join([mapping.get(g, g) for g in genes])
        
    # 5. Apply the conversion and save
    print("Replacing IDs with symbols...")
    df['intersections_symbols'] = df['intersections'].apply(replace_with_symbols)
    
    df.to_csv(output_csv, index=False)
    print(f"Success! Saved to {output_csv}")

# Run the function on your uploaded file
input_file = "gProfiler_athaliana_3-13-2026_2-26-35 PM__intersections.csv"
output_file = "gProfiler_athaliana_with_symbols.csv"

convert_gprofiler_ids_to_symbols(input_file, output_file)