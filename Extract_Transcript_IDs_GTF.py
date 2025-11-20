import os

def extract_transcript_id(file_path):
    transcript_ids = set()

    with open(file_path, 'r') as file:
        for line in file:
            if 'transcript_id' in line:
                try:
                    transcript_id = line.split('gene_id "')[1].split('"')[0]
                    transcript_ids.add(transcript_id)
                except IndexError:
                    continue  # skip malformed lines

    return transcript_ids


def main():
    # -----------------------------------
    # Set your input/output files here
    # -----------------------------------
    input_file = "F:\\Sorghum_MetaDEG\\Isoformswith\\Sbicolor_Merge.gtf"         # <-- change this
    output_file = "F:\\Sorghum_MetaDEG\\Alternative Splicing\\Figures\\IsoFormPergene\\SbiMerge_gene.txt"  # <-- change this
    # -----------------------------------

    transcript_ids = extract_transcript_id(input_file)

    with open(output_file, 'w') as out_file:
        for transcript_id in sorted(transcript_ids):
            out_file.write(f"{transcript_id}\n")

    print(f"Saved {len(transcript_ids)} transcript IDs to: {output_file}")


if __name__ == "__main__":
    main()
