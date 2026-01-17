import os
import subprocess

# --- CONFIGURATION ---
# Paths to your reference files
REF_TEXT = "/home/pgb-lab/Documents/Sorghum_circRNA/Drought/GenePread/Sbi_merged_circ.refFlat"
GENOME_FA = "/home/pgb-lab/Documents/Sorghum_circRNA/sbi_STAR_Index/sbi.fasta"

# Specify your input folder (where the .junction files are)
INPUT_DIR = "/home/pgb-lab/Documents/Sorghum_circRNA/Drought/circRNA_Mapped/STAR_Output/Junction_Reads/"

# Specify your output folder (where you want results saved)
OUTPUT_DIR = "/home/pgb-lab/Documents/Sorghum_circRNA/Drought/circRNA_CircExp/Results"
# ---------------------

def run_circ_pipeline():
    # Create output directory if it doesn't exist
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
        print "Created output directory: " + OUTPUT_DIR

    # Get list of junction files
    try:
        files = os.listdir(INPUT_DIR)
    except OSError:
        print "Error: Could not find input directory " + INPUT_DIR
        return

    junction_files = [f for f in files if f.endswith('_Chimeric.out.junction')]
    
    if not junction_files:
        print "No junction files found in " + INPUT_DIR
        return

    print "Found " + str(len(junction_files)) + " samples. Starting processing...\n"

    for junction_file in junction_files:
        # Extract ID (e.g., SRR17287606)
        sample_id = junction_file.split('_')[0]
        
        # Define full paths for input and output
        input_path = os.path.join(INPUT_DIR, junction_file)
        
        back_spliced_bed = os.path.join(OUTPUT_DIR, sample_id + "_back_spliced_junction.bed")
        parse_log = os.path.join(OUTPUT_DIR, sample_id + "_CIRCexplorer2_parse.log")
        output_txt = os.path.join(OUTPUT_DIR, sample_id + "_circularRNA_known.txt")
        annotate_log = os.path.join(OUTPUT_DIR, sample_id + "_CIRCexplorer2_annotate.log")
        
        print "--- Processing Sample: " + sample_id + " ---"

        # STEP 1: PARSE
        parse_cmd = "CIRCexplorer2 parse -t STAR -b " + back_spliced_bed + " " + input_path + " > " + parse_log + " 2>&1"
        
        # STEP 2: ANNOTATE 
        annotate_cmd = "CIRCexplorer2 annotate -r " + REF_TEXT + " -g " + GENOME_FA + " -b " + back_spliced_bed + " -o " + output_txt + " > " + annotate_log + " 2>&1"

        try:
            print "  Running 'parse'..."
            subprocess.check_call(parse_cmd, shell=True)
            
            print "  Running 'annotate'..."
            subprocess.check_call(annotate_cmd, shell=True)
            
            print "  [SUCCESS] Results saved to " + output_txt + "\n"
        except subprocess.CalledProcessError:
            print "  [ERROR] Sample " + sample_id + " failed. Check logs in output directory.\n"

    print "--- All tasks finished ---"

if __name__ == "__main__":
    run_circ_pipeline()