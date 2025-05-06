import argparse
from typing import List, Tuple, Set

def parse_gtf(gtf_file: str) -> List[Tuple[str, str]]:
    """Parse GTF file and extract unique gene_id-transcript_id pairs from transcript features"""
    gene_transcripts: List[Tuple[str, str]] = []
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            feature_type = parts[2]

            if feature_type == 'transcript':
                attributes = {}
                attr_str = parts[8]
                for attr in attr_str.split(';'):
                    attr = attr.strip()
                    if not attr:
                        continue
                    key, value = attr.split(' ', 1)
                    attributes[key] = value.strip('"')

                gene_id = attributes.get('gene_id')
                transcript_id = attributes.get('transcript_id')
                if gene_id and transcript_id:
                    gene_transcripts.append((gene_id, transcript_id))

    # Remove duplicates while preserving order
    seen: Set[Tuple[str, str]] = set()
    return [x for x in gene_transcripts if not (x in seen or seen.add(x))]

def read_gene_list(filename: str) -> Set[str]:
    """Read gene IDs from a file into a set"""
    with open(filename, 'r') as f:
        return {line.strip() for line in f}

def main():
    parser = argparse.ArgumentParser(description='Extract gene-transcript pairs based on gene list')
    parser.add_argument('--gtf', required=True, help='Input GTF file')
    parser.add_argument('--genes', required=True, help='File containing gene IDs to filter')
    parser.add_argument('--output', required=True, help='Output file name')
    args = parser.parse_args()

    # Parse GTF and get unique gene-transcript pairs
    gt_pairs = parse_gtf(args.gtf)

    #Read target genes
    target_genes = read_gene_list(args.genes)

    # Filter pairs
    filtered_pairs = [(g, t) for g, t in gt_pairs if g in target_genes]

    # Write output
    with open(args.output, 'w') as fout:
        fout.write("Gene_ID\tTranscript_ID\n")
        for gene_id, transcript_id in filtered_pairs:
            fout.write(f"{gene_id}\t{transcript_id}\n")

if __name__ == '__main__':
    main()
