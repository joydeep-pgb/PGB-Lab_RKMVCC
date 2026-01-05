#!/usr/bin/env python3

import statistics
from collections import defaultdict

# ===============================
# CONFIGURATION
# ===============================
GTF_FILE = "K:\\Sorghum_MetaDEG\\Isoform_Identification\\Sorghum_SUPPA\\Sbicolor_Merge.gtf"
# ===============================


def parse_attributes(attr_str):
    """Parse GTF attributes into a dictionary"""
    attrs = {}
    for field in attr_str.strip().split(";"):
        field = field.strip()
        if field:
            key, value = field.split(" ", 1)
            attrs[key] = value.replace('"', '')
    return attrs


def gtf_statistics(gtf_file):
    genes = set()
    transcripts = set()

    gene_to_transcripts = defaultdict(set)
    transcript_to_exons = defaultdict(list)
    transcript_lengths = {}
    exon_lengths = []

    strand_count = defaultdict(int)

    with open(gtf_file) as fh:
        for line in fh:
            if line.startswith("#"):
                continue

            cols = line.rstrip().split("\t")
            if len(cols) != 9:
                continue

            chrom, source, feature, start, end, score, strand, frame, attrs = cols
            start, end = int(start), int(end)

            attr_dict = parse_attributes(attrs)

            gene_id = attr_dict.get("gene_id")
            transcript_id = attr_dict.get("transcript_id")

            if feature == "transcript":
                genes.add(gene_id)
                transcripts.add(transcript_id)
                gene_to_transcripts[gene_id].add(transcript_id)
                transcript_lengths[transcript_id] = end - start + 1
                strand_count[strand] += 1

            elif feature == "exon":
                transcript_to_exons[transcript_id].append((start, end))
                exon_lengths.append(end - start + 1)

    transcripts_per_gene = [len(v) for v in gene_to_transcripts.values()]
    exons_per_transcript = [len(v) for v in transcript_to_exons.values()]

    print("\n===== GTF SUMMARY STATISTICS =====\n")

    print(f"Total genes                : {len(genes)}")
    print(f"Total transcripts          : {len(transcripts)}")
    print(f"Total exons                : {sum(exons_per_transcript)}\n")

    print("---- Transcripts per gene ----")
    print(f"Min                         : {min(transcripts_per_gene)}")
    print(f"Max                         : {max(transcripts_per_gene)}")
    print(f"Mean                        : {statistics.mean(transcripts_per_gene):.2f}")
    print(f"Median                      : {statistics.median(transcripts_per_gene)}\n")

    print("---- Exons per transcript ----")
    print(f"Min                         : {min(exons_per_transcript)}")
    print(f"Max                         : {max(exons_per_transcript)}")
    print(f"Mean                        : {statistics.mean(exons_per_transcript):.2f}")
    print(f"Median                      : {statistics.median(exons_per_transcript)}\n")

    print("---- Transcript length (bp) ----")
    print(f"Min                         : {min(transcript_lengths.values())}")
    print(f"Max                         : {max(transcript_lengths.values())}")
    print(f"Mean                        : {statistics.mean(transcript_lengths.values()):.2f}")
    print(f"Median                      : {statistics.median(transcript_lengths.values())}\n")

    print("---- Exon length (bp) ----")
    print(f"Min                         : {min(exon_lengths)}")
    print(f"Max                         : {max(exon_lengths)}")
    print(f"Mean                        : {statistics.mean(exon_lengths):.2f}")
    print(f"Median                      : {statistics.median(exon_lengths)}\n")

    print("---- Strand distribution ----")
    for strand in strand_count:
        print(f"Strand {strand}              : {strand_count[strand]}")

    print("\n=================================\n")


if __name__ == "__main__":
    gtf_statistics(GTF_FILE)
