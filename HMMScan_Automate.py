#!/usr/bin/env python3
import subprocess
import sys
import os
import pandas as pd

# ====================== USER CONFIGURATION ======================
HMM_PROFILE = "/path/to/your/profiles.hmm"   # <- EDIT THIS
FASTA_FILE = "/path/to/your/query.fasta"     # <- EDIT THIS
OUTPUT_BASE = "hmmscan_results"              # Output prefix
EVALUE = 1e-5                                 # Sequence E-value cutoff
DOMEVALUE = 1e-3                              # Domain E-value cutoff
CPU = 4                                       # Number of CPU cores
EXTRA_OPTIONS = ["--notextw"]                 # Optional flags
# ================================================================

def check_hmmer_installed():
    try:
        subprocess.run(["hmmscan", "-h"], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except (subprocess.CalledProcessError, FileNotFoundError):
        sys.exit("Error: HMMER not found. Install with conda (`conda install -c bioconda hmmer`)")

def hmmpress_profile(hmm_file):
    required_exts = ['.h3m', '.h3i', '.h3f', '.h3p']
    if not all(os.path.exists(hmm_file + ext) for ext in required_exts):
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

    print(f"[i] Running hmmscan:\n{' '.join(cmd)}\n")
    with open(outputs['stdout'], 'w') as outfile:
        subprocess.run(cmd, check=True, stdout=outfile)
    print(f"[✓] Scan complete. Output written to:")
    for key, path in outputs.items():
        print(f"   - {key}: {path}")

    return outputs

def parse_tblout(tbl_file, out_csv):
    """Parse --tblout file (sequence level hits) into a DataFrame."""
    print(f"[i] Parsing tblout: {tbl_file}")
    records = []
    with open(tbl_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split()
            if len(fields) >= 18:
                records.append({
                    "target_name": fields[0],
                    "target_accession": fields[1],
                    "query_name": fields[2],
                    "query_accession": fields[3],
                    "E-value": float(fields[4]),
                    "score": float(fields[5]),
                    "bias": float(fields[6]),
                    "exp": float(fields[9]),
                    "description": " ".join(fields[18:])
                })
    df = pd.DataFrame(records)
    df.to_csv(out_csv, index=False)
    print(f"[✓] Parsed tblout saved to: {out_csv}")
    return df

def parse_domtblout(domtbl_file, out_csv):
    """Parse --domtblout file (domain hits) into a DataFrame."""
    print(f"[i] Parsing domtblout: {domtbl_file}")
    records = []
    with open(domtbl_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split()
            if len(fields) >= 23:
                records.append({
                    "target_name": fields[0],
                    "target_accession": fields[1],
                    "query_name": fields[3],
                    "E-value": float(fields[6]),
                    "score": float(fields[7]),
                    "bias": float(fields[8]),
                    "domain_coord": f"{fields[17]}-{fields[18]}",
                    "description": " ".join(fields[22:])
                })
    df = pd.DataFrame(records)
    df.to_csv(out_csv, index=False)
    print(f"[✓] Parsed domtblout saved to: {out_csv}")
    return df

def main():
    if not os.path.exists(HMM_PROFILE):
        sys.exit(f"HMM profile not found: {HMM_PROFILE}")
    if not os.path.exists(FASTA_FILE):
        sys.exit(f"FASTA file not found: {FASTA_FILE}")

    check_hmmer_installed()
    hmmpress_profile(HMM_PROFILE)
    output_files = run_hmmscan()

    # Parsing outputs
    parse_tblout(output_files['tblout'], OUTPUT_BASE + ".tblout.parsed.csv")
    parse_domtblout(output_files['domtblout'], OUTPUT_BASE + ".domtblout.parsed.csv")

if __name__ == "__main__":
    main()
