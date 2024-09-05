import pandas as pd
import matplotlib.pyplot as plt

def calculate_gene_density_profile(df, bin_size, window_size):
    max_position = df['End'].max()
    bins = range(0, max_position, bin_size)
    gene_density = []

    for start in bins:
        end = start + window_size
        genes_in_window = ((df['Start'] >= start) & (df['End'] <= end)).sum()
        gene_density.append(genes_in_window)

    return gene_density

# Define the file paths and parameters
input_file = r'D:\Bioinformatics_Data\Sorghum_lncRNA\Classification\Gene density\Others\other_lncRNA_cls.txt'
output_prefix = 'gene_density_profile'
bin_size = 100000
window_size = 1000000

# Read the TSV file into a DataFrame (skip header)
df = pd.read_csv(input_file, sep='\t', skiprows=1, header=None, names=['GeneID', 'Chromosome', 'Start', 'End'])

# Get unique chromosomes
unique_chromosomes = df['Chromosome'].unique()

# Calculate and save gene density profile for each chromosome
for chromosome in unique_chromosomes:
    chromosome_df = df[df['Chromosome'] == chromosome]
    gene_density = calculate_gene_density_profile(chromosome_df, bin_size, window_size)

    # Save the gene density profile to a TSV file
    output_file = f'{output_prefix}_{chromosome}.tsv'
    bin_positions = range(0, len(gene_density) * bin_size, bin_size)
    df_output = pd.DataFrame({'Bin Position': bin_positions, 'Gene Density': gene_density})
    df_output.to_csv(output_file, sep='\t', index=False)

    # Plot the gene density profile
    plt.plot(bin_positions, gene_density)
    plt.xlabel('Genomic Position')
    plt.ylabel('Gene Density')
    plt.title(f'Gene Density Profile - Chromosome {chromosome}')
    plt.savefig(f'gene_density_profile_chrom_{chromosome}.png')
    plt.clf()  # Clear the plot for the next chromosome

# Show a confirmation message
print("Gene density profiles saved successfully.")
