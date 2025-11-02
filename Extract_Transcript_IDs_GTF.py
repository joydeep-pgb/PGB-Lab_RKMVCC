import os

def extract_transcript_id(file_path):
    transcript_ids = set()

    with open(file_path, 'r') as file:
        for line in file:
            if 'transcript_id' in line:
                transcript_id = line.split('transcript_id "')[1].split('"')[0]
                transcript_ids.add(transcript_id)

    return transcript_ids

def main():
    # === Set file paths here ===
    input_gtf = "/run/media/joydeep/One_HDD/Sorghum_MetaDEG/WGCNA/Sbicolor_DAS_Isoforms.gtf"
    output_file = "/run/media/joydeep/One_HDD/Sorghum_MetaDEG/WGCNA/DAS_IsoformIDs.txt"
    # ============================

    if not os.path.exists(input_gtf):
        print(f"Error: Input file not found -> {input_gtf}")
        return

    transcript_ids = extract_transcript_id(input_gtf)

    with open(output_file, 'w') as out_file:
        for transcript_id in sorted(transcript_ids):
            out_file.write(f"{transcript_id}\n")

    print(f"✅ Extracted {len(transcript_ids)} transcript IDs saved to: {output_file}")

if __name__ == "__main__":
    main()
