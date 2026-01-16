#!/usr/bin/env python3

import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import os

# ==============================
# 🔧 USER SETTINGS
# ==============================
BAM_DIR = "/run/media/joydeep/One_HDD/Sorghum_MetaDEG/Salt/Mapped/BAM_Files/"   # ← CHANGE THIS
MAX_WORKERS = min(8, os.cpu_count())           # adjust if needed
# ==============================


def index_bam(bam):
    """Index a single BAM file using samtools"""
    bai = bam.with_suffix(".bam.bai")

    if bai.exists():
        return f"[SKIP] {bam.name}"

    try:
        subprocess.run(
            ["samtools", "index", str(bam)],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )
        return f"[DONE] {bam.name}"
    except subprocess.CalledProcessError:
        return f"[ERROR] {bam.name}"


def main():
    bam_dir = Path(BAM_DIR)

    if not bam_dir.is_dir():
        raise SystemExit(f"[ERROR] Directory not found: {bam_dir}")

    bam_files = sorted(bam_dir.glob("*.bam"))

    if not bam_files:
        print("[WARNING] No BAM files found.")
        return

    print(f"[INFO] Found {len(bam_files)} BAM files")
    print(f"[INFO] Using {MAX_WORKERS} threads\n")

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = [executor.submit(index_bam, bam) for bam in bam_files]

        for future in as_completed(futures):
            print(future.result())

    print("\n[DONE] Multithreaded BAM indexing completed.")


if __name__ == "__main__":
    main()
