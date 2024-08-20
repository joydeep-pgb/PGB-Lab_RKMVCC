import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler

# Function to read the gene expression data
def read_data(file_path):
    # Assuming the data is tab-separated
    data = pd.read_csv(file_path, sep="\t", index_col=0)
    return data

# Function to find the optimal number of K using the elbow method
def find_optimal_k(data, max_k=10):
    sse = []
    k_range = range(1, max_k + 1)
    
    for k in k_range:
        kmeans = KMeans(n_clusters=k, n_init=10, random_state=42)
        kmeans.fit(data)
        sse.append(kmeans.inertia_)
    
    plt.figure(figsize=(8, 5))
    plt.plot(k_range, sse, 'bx-')
    plt.xlabel('Number of clusters K')
    plt.ylabel('Sum of squared distances (SSE)')
    plt.title('Elbow Method For Optimal K')
    plt.show()

# Main function to execute the analysis
def main(file_path):
    # Load the data
    data = read_data(file_path)
    
    # Standardize the data (important for gene expression data)
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(data)
    
    # Find the optimal number of clusters
    find_optimal_k(scaled_data)

if __name__ == "__main__":
    # Replace 'your_gene_expression_data.txt' with the path to your file
    file_path = 'K-Means.txt'
    main(file_path)
