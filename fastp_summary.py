import os
import pandas as pd
from bs4 import BeautifulSoup

# Function to extract the total bases and total reads from an HTML file
def extract_data_from_html(file_path):
    with open(file_path, 'r') as file:
        soup = BeautifulSoup(file, 'html.parser')
        
        # Extracting the total reads and total bases before filtering
        before_filtering_section = soup.find(id='before_filtering_summary')
        total_reads_before = before_filtering_section.find(text='total reads:').find_next('td').text.strip()
        total_bases_before = before_filtering_section.find(text='total bases:').find_next('td').text.strip()
        
        # Extracting the total reads and total bases after filtering
        after_filtering_section = soup.find(id='after_filtering_summary')
        total_reads_after = after_filtering_section.find(text='total reads:').find_next('td').text.strip()
        total_bases_after = after_filtering_section.find(text='total bases:').find_next('td').text.strip()
        
        return total_reads_before, total_bases_before, total_reads_after, total_bases_after

# Function to process all HTML files in a directory
def process_html_files(directory_path, output_file):
    data = []
    
    # Iterate through each file in the directory
    for filename in os.listdir(directory_path):
        if filename.endswith('.html'):
            file_path = os.path.join(directory_path, filename)
            
            # Extract data from the HTML file
            total_reads_before, total_bases_before, total_reads_after, total_bases_after = extract_data_from_html(file_path)
            
            # Append the data to the list
            data.append([filename, total_reads_before, total_bases_before, total_reads_after, total_bases_after])
    
    # Create a DataFrame
    df = pd.DataFrame(data, columns=['Sample Name', 'Total Reads Before', 'Total Bases Before', 'Total Reads After', 'Total Bases After'])
    
    # Save the DataFrame to a CSV file
    df.to_csv(output_file, sep='\t', index=False)
    
    # Optional: Print confirmation
    print(f"Data has been saved to '{output_file}'.")

# Example usage of the function
directory_path = 'fastp_LF/'  # Replace with the path to your HTML files
output_file = 'Table3.txt'  # Specify the output file name

process_html_files(directory_path, output_file)
