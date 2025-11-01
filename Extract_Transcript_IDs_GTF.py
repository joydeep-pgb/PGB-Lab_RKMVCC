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
    # ✅ Set your input and output file paths here
    input_gtf = "/media/pgb-lab/One_HDD/Sorghum_MetaDEG/Isoformswith/Drought/Sbicolor_Drought_DAS.gtf"
    output_file = "/media/pgb-lab/One_HDD/Sorghum_MetaDEG/Isoformswith/Drought/Transcript_IDs_Sbicolor_Drought_DAS.txt"

    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    transcript_ids = extract_transcript_id(input_gtf)

    with open(output_file, 'w') as out_file:
        for transcript_id in sorted(transcript_ids):
            out_file.write(f"{transcript_id}\n")

    print(f"✅ Extracted {len(transcript_ids)} transcript IDs.")
    print(f"💾 Saved to: {output_file}")


if __name__ == "__main__":
    main()
