#!/usr/bin/env python3

def filter_gtf(gtf_file, id_file, output_file):
    """
    Extracts lines from a GTF file where gene_id exactly matches
    any ID listed in the provided file.

    Parameters:
    -----------
    gtf_file : str
        Path to input GTF file.
    id_file : str
        Text file containing one gene_id per line (e.g., MSTRG.13).
    output_file : str
        Path to write the filtered GTF.
    """

    # Load target gene IDs
    with open(id_file, 'r') as f:
        target_ids = set(line.strip() for line in f if line.strip())

    print(f"Loaded {len(target_ids)} target gene IDs from {id_file}")

    with open(gtf_file, 'r') as infile, open(output_file, 'w') as outfile:
        kept = 0
        total = 0

        for line in infile:
            total += 1

            # Keep header lines
            if line.startswith('#'):
                outfile.write(line)
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            attributes = fields[8]

            # Extract gene_id
            gene_id = None
            for attr in attributes.split(';'):
                attr = attr.strip()
                if attr.startswith('gene_id'):
                    gene_id = attr.split('"')[1]
                    break

            # Write matching lines
            if gene_id in target_ids:
                outfile.write(line)
                kept += 1

        print(f"Processed {total} lines.")
        print(f"Kept {kept} lines matching provided gene IDs.")
        print(f"Filtered GTF saved to: {output_file}")


# ------------------------
# 🔧 SET YOUR FILE PATHS HERE
# ------------------------

gtf_file = "/run/media/joydeep/One_HDD/Mango_AS_rMATS/rMATS/mango_AS_meta.gtf"
id_file = "/run/media/joydeep/One_HDD/Mango_AS_rMATS/rMATS/rMATS DATA/Final Events/Total_DAS_Events.txt"
output_file = "/run/media/joydeep/One_HDD/Mango_AS_rMATS/rMATS/DAS_filtered.gtf"

# ------------------------
# Run filtering
# ------------------------
if __name__ == "__main__":
    filter_gtf(gtf_file, id_file, output_file)
