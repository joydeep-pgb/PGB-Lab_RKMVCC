#!/bin/bash

# Main folder where GTF files are located
MAIN_FOLDER="/home/pgb-lab/Documents/Sorghum_MetaDEG/Drought/Assembled/"  # Change this to your folder path
WORKING_DIR="$MAIN_FOLDER/Working_directory"

# Create Working_directory if not exists
mkdir -p "$WORKING_DIR"

# Create Counts folder inside Working_directory
mkdir -p "$WORKING_DIR/Counts"

# Loop through all GTF files in the main folder
for gtf_file in "$MAIN_FOLDER"/*.gtf; do
    # Get the filename without the path
    filename=$(basename "$gtf_file")
    
    # Extract SRA name without extension
    sra_name="${filename%.gtf}"

    # Copy GTF file to Working_directory
    cp "$gtf_file" "$WORKING_DIR"

    # Create folder inside Working_directory with SRA name
    mkdir -p "$WORKING_DIR/$sra_name"

    # Move the copied GTF file to the corresponding folder
    mv "$WORKING_DIR/$filename" "$WORKING_DIR/$sra_name/"

    # Change to Working_directory
    cd "$WORKING_DIR" || exit

    # Run prepDE.py command with appropriate filenames
    prepDE.py -l 75 -g "${sra_name}_gene.csv" -t "${sra_name}_transcript.csv"

    # Move output CSV files to Counts folder inside Working_directory
    mv "${sra_name}_gene.csv" "${sra_name}_transcript.csv" "$WORKING_DIR/Counts/"

    # Delete SRA folder after processing
    rm -r "$WORKING_DIR/$sra_name"

    # Return to the main directory
    cd "$MAIN_FOLDER" || exit

done

echo "All GTF files have been processed, and counts have been moved to the Counts folder."
