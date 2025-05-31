#!/usr/bin/env python3
import subprocess
import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from dna_features_viewer import GraphicFeature, GraphicRecord

# ====================== USER CONFIGURATION ======================
HMM_PROFILE = "/path/to/Pfam-A.hmm"              # <- EDIT THIS
FASTA_FILE = "/path/to/proteome.fasta"           # <- EDIT THIS
OUTPUT_BASE = "abc_hmmscan"
PFAM_IDS = {"PF00005", "PF00664", "PF01061"}     # <- Add more Pfam domains
EVALUE = 1e-5
DOMEVALUE = 1e-3
CPU = 4
EXTRA_OPTIONS = ["--notextw"]
PLOT_DOMAIN_COUNTS = True
PLOT_DOMAIN_ARCHITECTURE = True
# ================================================================

def check_hmmer_installed():
    try:
        subprocess.run(["hmmscan", "-h"], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except Exception:
        sys.exit("Error: HMMER not found. Install with: conda install -c bioconda hmmer")

def hmmpress_profile(hmm_file):
    exts = ['.h3m', '.h3i', '.h3f', '.h3p']
    if not all(os.path.exists(hmm_file + ext) for ext in exts):
        print(f"[i] Pressing HMM profile: {hmm_file}")
        subprocess.run(["hmmpress", hmm_file], check=True)
    else:
        print("[✓] HMM profile already pressed.")

def run_hmmscan():
    outputs = {
        'domtblout': f"{OUTPUT_BASE}.domtblout",
        'tblout': f"{OUTPUT_BASE}.tblout",
        'stdout': f"{OUTPUT_BASE}.full.txt"
    }
    cmd = [
        "hmmscan",
        "--cpu", str(CPU),
        "-E", str(EVALUE),
        "--domE", str(DOMEVALUE),
        "--tblout", outputs['tblout'],
        "--domtblout", outputs['domtblout']
    ] + EXTRA_OPTIONS + [HMM_PROFILE, FASTA_FILE]
    print(f"[i] Running hmmscan...")
    with open(outputs['stdout'], 'w') as out:
        subprocess.run(cmd, check=True, stdout=out)
    return outputs

def parse_and_filter(domtbl_file, filtered_csv, id_txt):
    records = []
    protein_ids = set()

    with open(domtbl_file) as f:
        for line in f:
            if line.startswith("#"): continue
            fields = line.strip().split()
            if len(fields) >= 23:
                pfam_id = fields[1]
                if pfam_id in PFAM_IDS:
                    pid = fields[0]
                    records.append({
                        "protein_id": pid,
                        "pfam_id": pfam_id,
                        "query": fields[3],
                        "evalue": float(fields[6]),
                        "score": float(fields[7]),
                        "bias": float(fields[8]),
                        "start": int(fields[17]),
                        "end": int(fields[18]),
                        "description": " ".join(fields[22:])
                    })
                    protein_ids.add(pid)

    df = pd.DataFrame(records)
    df.to_csv(filtered_csv, index=False)
    with open(id_txt, "w") as f:
        for pid in sorted(protein_ids):
            f.write(pid + "\n")
    print(f"[✓] {len(df)} filtered hits saved to: {filtered_csv}")
    return df

def plot_domain_counts(df, outfile):
    plt.figure(figsize=(8, 5))
    sns.countplot(data=df, x='pfam_id', order=df['pfam_id'].value_counts().index)
    plt.ylabel("Domain Occurrences")
    plt.title("Pfam Domain Frequency")
    plt.tight_layout()
    plt.savefig(outfile)
    print(f"[✓] Domain frequency plot saved to: {outfile}")

def plot_clean_architectures(df, outfile, max_proteins=20):
    proteins = df['protein_id'].unique()[:max_proteins]
    records = []

    for pid in proteins:
        features = []
        p_df = df[df['protein_id'] == pid]
        for _, row in p_df.iterrows():
            color = "#4CAF50" if "membrane" in row["description"].lower() else "#3F51B5"
            features.append(GraphicFeature(
                start=row['start'],
                end=row['end'],
                strand=+1,
                color=color,
                label=row['pfam_id']
            ))
        length = max(p_df['end']) + 50
        record = GraphicRecord(sequence_length=length, features=features)
        records.append((pid, record))

    fig, axs = plt.subplots(len(records), 1, figsize=(10, 1.5 * len(records)))
    if len(records) == 1:
        axs = [axs]

    for ax, (pid, record) in zip(axs, records):
        record.plot(ax=ax)
        ax.set_title(pid, loc='left', fontsize=9)

    plt.tight_layout()
    plt.savefig(outfile)
    print(f"[✓] Clean domain layout saved to: {outfile}")

def main():
    if not os.path.exists(HMM_PROFILE): sys.exit(f"HMM file not found: {HMM_PROFILE}")
    if not os.path.exists(FASTA_FILE): sys.exit(f"FASTA file not found: {FASTA_FILE}")

    check_hmmer_installed()
    hmmpress_profile(HMM_PROFILE)
    files = run_hmmscan()

    filtered_csv = OUTPUT_BASE + ".ABC_filtered.csv"
    id_txt = OUTPUT_BASE + ".ABC_ids.txt"
    df = parse_and_filter(files['domtblout'], filtered_csv, id_txt)

    if PLOT_DOMAIN_COUNTS:
        plot_domain_counts(df, OUTPUT_BASE + "_domain_counts.png")
    if PLOT_DOMAIN_ARCHITECTURE:
        plot_clean_architectures(df, OUTPUT_BASE + "_clean_architecture.png")

if __name__ == "__main__":
    main()
