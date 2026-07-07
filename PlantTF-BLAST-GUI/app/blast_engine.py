"""
blast_engine.py
----------------
Non-GUI logic for the PlantTF-BLAST application:

  * checking that NCBI BLAST+ is installed
  * detecting protein vs. nucleotide FASTA
  * building (and caching) a BLAST protein database from a user-supplied
    "custom TF database" FASTA file whose headers follow the PlantTFDB
    convention:

        >GeneID Species name|TF_Family|Free text description

  * running blastp / blastx in chunks (so the GUI can show progress and
    the user can cancel a long run)
  * parsing raw BLAST tabular output into structured TFHit records
  * collapsing raw HSPs into one best record per query, PlantTFDB-style

Nothing in this file touches PyQt; it can be imported and unit-tested
on its own.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable, List, Optional, Sequence

# --------------------------------------------------------------------------
# Data structures
# --------------------------------------------------------------------------


@dataclass
class DBEntry:
    """One parsed header line from the custom TF database FASTA."""

    seq_id: str
    species: str
    family: str
    description: str
    raw_header: str


@dataclass
class TFHit:
    """A single BLAST HSP, annotated with TF family/description info."""

    query_id: str
    hit_id: str
    family: str
    species: str
    description: str
    evalue: float
    bitscore: float
    pident: float
    qcov: float
    length: int


@dataclass
class QueryResult:
    """Aggregated, PlantTFDB-style result for one query sequence."""

    query_id: str
    matched: bool
    family: str = "-"
    hit_id: str = "-"
    evalue: Optional[float] = None
    bitscore: Optional[float] = None
    pident: Optional[float] = None
    qcov: Optional[float] = None
    description: str = ""
    species: str = ""
    all_hits: List[TFHit] = field(default_factory=list)


class BlastError(RuntimeError):
    """Raised for any recoverable BLAST/pipeline failure."""


class BlastCancelled(RuntimeError):
    """Raised internally when the user cancels a running job."""


# --------------------------------------------------------------------------
# Dependency / environment checks
# --------------------------------------------------------------------------

REQUIRED_TOOLS = ("makeblastdb", "blastp", "blastx")


def find_blast_tools() -> dict:
    """Return {tool_name: path_or_None} for the BLAST+ executables we need."""
    return {tool: shutil.which(tool) for tool in REQUIRED_TOOLS}


def blast_is_installed() -> bool:
    return all(find_blast_tools().values())


def install_hint() -> str:
    """A short, platform-aware suggestion for installing NCBI BLAST+."""
    return (
        "NCBI BLAST+ was not found on your PATH.\n\n"
        "Install it with one of the following, then restart this app:\n\n"
        "  Fedora / RHEL:   sudo dnf install ncbi-blast+\n"
        "  Debian / Ubuntu: sudo apt install ncbi-blast+\n"
        "  Conda:           conda install -c bioconda blast\n"
        "  macOS (brew):    brew install blast\n\n"
        "Or download binaries from:\n"
        "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/"
    )


# --------------------------------------------------------------------------
# FASTA helpers
# --------------------------------------------------------------------------

_AMINO_ONLY = set("ACDEFGHIKLMNPQRSTVWY")
_NUCLEOTIDE = set("ACGTUN")


def detect_sequence_type(fasta_path: str, sample_seqs: int = 20) -> str:
    """Return 'protein' or 'nucleotide' by sampling the first few sequences."""
    seen = 0
    letters = []
    with open(fasta_path, "r", errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                seen += 1
                if seen > sample_seqs:
                    break
                continue
            letters.append(line.strip().upper())
    joined = "".join(letters)[:5000]
    if not joined:
        return "protein"
    non_nuc = sum(1 for c in joined if c not in _NUCLEOTIDE and c.isalpha())
    ratio_non_nuc = non_nuc / max(1, len(joined))
    # if more than ~8% of characters are outside ACGTUN it is protein
    return "protein" if ratio_non_nuc > 0.08 else "nucleotide"


def count_fasta_sequences(fasta_path: str) -> int:
    count = 0
    with open(fasta_path, "r", errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                count += 1
    return count


def parse_db_header(header_line: str) -> DBEntry:
    """
    Parse one '>' header from the custom TF database FASTA.

    Expected PlantTFDB-style convention:
        >GeneID Species name|Family|Description

    Falls back gracefully if the pipe-delimited annotation is missing.
    """
    raw = header_line[1:].strip() if header_line.startswith(">") else header_line.strip()
    if " " in raw:
        seq_id, rest = raw.split(" ", 1)
    else:
        seq_id, rest = raw, ""

    parts = rest.split("|")
    if len(parts) >= 3:
        species = parts[0].strip()
        family = parts[1].strip()
        description = "|".join(parts[2:]).strip()
    elif len(parts) == 2:
        species = ""
        family = parts[0].strip()
        description = parts[1].strip()
    else:
        species = ""
        family = "Unknown"
        description = rest.strip()

    return DBEntry(
        seq_id=seq_id,
        species=species,
        family=family or "Unknown",
        description=description,
        raw_header=raw,
    )


def parse_all_db_headers(fasta_path: str) -> dict:
    """Return {seq_id: DBEntry} for every sequence in the TF database FASTA."""
    entries = {}
    with open(fasta_path, "r", errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                entry = parse_db_header(line)
                entries[entry.seq_id] = entry
    return entries


def split_fasta(fasta_path: str, chunk_size: int, out_dir: str) -> List[str]:
    """Split a FASTA file into chunk files of `chunk_size` sequences each.
    Returns list of chunk file paths (all written into out_dir)."""
    chunk_paths = []
    buf: List[str] = []
    n_in_chunk = 0
    chunk_idx = 0

    def flush():
        nonlocal buf, n_in_chunk, chunk_idx
        if not buf:
            return
        chunk_idx += 1
        path = os.path.join(out_dir, f"chunk_{chunk_idx:05d}.fasta")
        with open(path, "w") as out:
            out.writelines(buf)
        chunk_paths.append(path)
        buf = []
        n_in_chunk = 0

    with open(fasta_path, "r", errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                if n_in_chunk >= chunk_size:
                    flush()
                n_in_chunk += 1
            buf.append(line)
    flush()
    return chunk_paths


# --------------------------------------------------------------------------
# BLAST database caching
# --------------------------------------------------------------------------


def get_cache_root() -> Path:
    root = Path(tempfile.gettempdir()) / "planttf_blast_gui_cache"
    root.mkdir(parents=True, exist_ok=True)
    return root


def _fasta_signature(fasta_path: str) -> str:
    st = os.stat(fasta_path)
    key = f"{os.path.abspath(fasta_path)}|{st.st_size}|{st.st_mtime_ns}"
    return hashlib.sha1(key.encode("utf-8")).hexdigest()[:16]


def ensure_blast_db(
    fasta_path: str,
    log: Optional[Callable[[str], None]] = None,
    force_rebuild: bool = False,
) -> str:
    """
    Build (or reuse a cached) protein BLAST database for `fasta_path`.
    Returns the db prefix path suitable for `-db`.
    """
    log = log or (lambda msg: None)
    sig = _fasta_signature(fasta_path)
    db_dir = get_cache_root() / sig
    db_prefix = str(db_dir / "db")
    marker = db_dir / "signature.json"

    if not force_rebuild and marker.exists():
        try:
            info = json.loads(marker.read_text())
            if info.get("source") == os.path.abspath(fasta_path) and (
                db_dir / "db.phr"
            ).exists():
                log(f"Using cached BLAST database ({db_dir.name}).")
                return db_prefix
        except Exception:
            pass

    if db_dir.exists():
        shutil.rmtree(db_dir, ignore_errors=True)
    db_dir.mkdir(parents=True, exist_ok=True)

    log("Building BLAST database from the custom TF FASTA file...")
    cmd = [
        "makeblastdb",
        "-in",
        fasta_path,
        "-dbtype",
        "prot",
        "-out",
        db_prefix,
        "-parse_seqids",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        # -parse_seqids can fail on duplicate/odd IDs; retry without it
        log("Retrying database build without -parse_seqids ...")
        cmd = ["makeblastdb", "-in", fasta_path, "-dbtype", "prot", "-out", db_prefix]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise BlastError(f"makeblastdb failed:\n{proc.stderr or proc.stdout}")

    marker.write_text(json.dumps({"source": os.path.abspath(fasta_path)}))
    log("BLAST database ready.")
    return db_prefix


# --------------------------------------------------------------------------
# Running BLAST
# --------------------------------------------------------------------------

OUTFMT_FIELDS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "evalue",
    "bitscore",
    "qcovs",
    "stitle",
]
OUTFMT = "6 " + " ".join(OUTFMT_FIELDS)


def _parse_tabular_line(line: str) -> Optional[dict]:
    parts = line.rstrip("\n").split("\t")
    if len(parts) < len(OUTFMT_FIELDS):
        return None
    row = dict(zip(OUTFMT_FIELDS, parts))
    return row


def run_blast_chunked(
    query_path: str,
    db_prefix: str,
    program: str = "blastp",
    evalue: float = 1e-5,
    max_target_seqs: int = 5,
    threads: int = 4,
    chunk_size: int = 200,
    progress_cb: Optional[Callable[[int, int, str], None]] = None,
    should_stop: Optional[Callable[[], bool]] = None,
    log: Optional[Callable[[str], None]] = None,
) -> List[dict]:
    """
    Run BLAST on `query_path` against `db_prefix`, splitting the query into
    chunks so we can report progress and honour cancellation.

    Returns a list of raw parsed tabular rows (dicts).
    """
    log = log or (lambda msg: None)
    progress_cb = progress_cb or (lambda done, total, msg: None)
    should_stop = should_stop or (lambda: False)

    total_queries = count_fasta_sequences(query_path)
    if total_queries == 0:
        raise BlastError("The query FASTA file appears to contain no sequences.")

    work_dir = tempfile.mkdtemp(prefix="planttf_blast_run_")
    try:
        chunks = split_fasta(query_path, chunk_size, work_dir)
        log(f"Query split into {len(chunks)} chunk(s) ({total_queries} sequences total).")

        all_rows: List[dict] = []
        done_seqs = 0

        for i, chunk_path in enumerate(chunks, start=1):
            if should_stop():
                raise BlastCancelled("BLAST run was cancelled by the user.")

            n_in_chunk = count_fasta_sequences(chunk_path)
            cmd = [
                program,
                "-query",
                chunk_path,
                "-db",
                db_prefix,
                "-evalue",
                str(evalue),
                "-max_target_seqs",
                str(max(1, max_target_seqs)),
                "-num_threads",
                str(max(1, threads)),
                "-outfmt",
                OUTFMT,
            ]
            proc = subprocess.Popen(
                cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
            )
            stdout, stderr = proc.communicate()

            if proc.returncode != 0:
                raise BlastError(f"{program} failed on chunk {i}:\n{stderr}")

            if stderr:
                for line in stderr.splitlines():
                    line = line.strip()
                    if line:
                        log(f"[{program}] {line}")

            for line in stdout.splitlines():
                row = _parse_tabular_line(line)
                if row:
                    all_rows.append(row)

            done_seqs += n_in_chunk
            progress_cb(
                min(done_seqs, total_queries),
                total_queries,
                f"Searched {min(done_seqs, total_queries)}/{total_queries} query sequences...",
            )

        return all_rows
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)


# --------------------------------------------------------------------------
# Assembling PlantTFDB-style results
# --------------------------------------------------------------------------


def build_query_results(
    query_path: str,
    raw_rows: Sequence[dict],
    db_entries: dict,
    max_hits_per_query: int = 5,
) -> List[QueryResult]:
    """
    Collapse raw BLAST tabular rows into one QueryResult per query sequence
    (present in query_path, even if it had zero hits), keeping up to
    `max_hits_per_query` distinct best subject hits per query, ranked by
    e-value / bitscore, PlantTFDB-style (best hit = predicted TF family).
    """
    # group rows by query, de-duplicating repeated HSPs against the same
    # subject (keep only the first / best-scoring HSP per subject).
    per_query_best_hsp: dict = {}
    for row in raw_rows:
        qid = row["qseqid"]
        sid = row["sseqid"]
        try:
            evalue = float(row["evalue"])
            bitscore = float(row["bitscore"])
            pident = float(row["pident"])
            qcov = float(row["qcovs"]) if row["qcovs"] not in ("", "N/A") else 0.0
            length = int(row["length"])
        except ValueError:
            continue

        bucket = per_query_best_hsp.setdefault(qid, {})
        existing = bucket.get(sid)
        if existing is None or bitscore > existing.bitscore:
            entry = db_entries.get(sid)
            family = entry.family if entry else "Unknown"
            species = entry.species if entry else ""
            description = entry.description if entry else row.get("stitle", "")
            bucket[sid] = TFHit(
                query_id=qid,
                hit_id=sid,
                family=family,
                species=species,
                description=description,
                evalue=evalue,
                bitscore=bitscore,
                pident=pident,
                qcov=qcov,
                length=length,
            )

    # enumerate every query in the original FASTA so "no hit" queries appear too
    all_query_ids = []
    with open(query_path, "r", errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                qid = line[1:].strip().split()[0]
                all_query_ids.append(qid)

    results: List[QueryResult] = []
    for qid in all_query_ids:
        hits_dict = per_query_best_hsp.get(qid, {})
        hits = sorted(hits_dict.values(), key=lambda h: (h.evalue, -h.bitscore))
        top_hits = hits[:max_hits_per_query]

        if not top_hits:
            results.append(QueryResult(query_id=qid, matched=False))
            continue

        best = top_hits[0]
        results.append(
            QueryResult(
                query_id=qid,
                matched=True,
                family=best.family,
                hit_id=best.hit_id,
                evalue=best.evalue,
                bitscore=best.bitscore,
                pident=best.pident,
                qcov=best.qcov,
                description=best.description,
                species=best.species,
                all_hits=top_hits,
            )
        )
    return results


def family_distribution(results: Iterable[QueryResult]) -> "list[tuple[str, int]]":
    """Return (family, count) pairs, sorted by count descending, for matched queries."""
    counts: dict = {}
    for r in results:
        if r.matched:
            counts[r.family] = counts.get(r.family, 0) + 1
    return sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
