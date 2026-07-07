# PlantTF-BLAST GUI

A PyQt6 desktop application that lets you run a **custom BLAST search** against
your own species-specific transcription factor (TF) database and predicts TF
identity for query protein sequences, reporting results in the same style as
[PlantTFDB](http://planttfdb.gao-lab.org/)'s BLAST tool:

```
Query_ID              TF_Family   Best_Hit_ID       E-value     Description
Phvul.001G031300       bHLH        AT1G10120.1       3e-84       bHLH family protein
Phvul.001G035600       GATA        AT4G26150.1       4e-22       cytokinin-responsive gata factor 1
Phvul.001G037000       HSF         AT2G26150.1       1e-115      heat shock transcription factor A2
Phvul.001G042200       WRKY        AT1G80840.1       5e-76       WRKY DNA-binding protein 40
```

## How it works

1. You provide a **custom TF database FASTA** whose headers follow the
   PlantTFDB convention:

   ```
   >GeneID Species name|TF_Family|Description
   ```

   For example (rice, from PlantTFDB):

   ```
   >LOC_Os02g58440.1 Oryza sativa subsp. japonica|C3H|C3H family protein
   ```

   This is exactly the format PlantTFDB uses for its own downloadable
   per-species TF protein FASTA files, so you can use any PlantTFDB species
   download directly as your custom database.

2. You provide a **query FASTA** (protein, or nucleotide for blastx) — e.g.
   your differentially expressed gene set.

3. The app builds (and caches) a local BLAST protein database from your TF
   FASTA, then runs `blastp` (or `blastx` for nucleotide queries) in chunks
   in a background thread so the GUI stays responsive and the run can be
   cancelled.

4. For each query sequence, the best-scoring hit's TF family (parsed straight
   from the database header) is reported as the **predicted TF family**,
   along with the hit ID, e-value, %identity, query coverage, and
   description — mirroring PlantTFDB's own BLAST result table.

> **Note:** This is a homology-based (BLAST) classification, not a Pfam/HMM
> domain-based classification like PlantTFDB's internal pipeline. Results
> are only as good as the family annotations already present in your
> database headers.

## Requirements

- Python 3.9+
- [PyQt6](https://pypi.org/project/PyQt6/)
- **NCBI BLAST+** (`makeblastdb`, `blastp`, `blastx`) available on your `PATH`
- (optional) `openpyxl`, only if you want to export results as `.xlsx`

### Installing NCBI BLAST+

| OS | Command |
|---|---|
| Fedora / RHEL | `sudo dnf install ncbi-blast+` |
| Debian / Ubuntu | `sudo apt install ncbi-blast+` |
| Conda | `conda install -c bioconda blast` |
| macOS (Homebrew) | `brew install blast` |

Or download prebuilt binaries from the
[NCBI BLAST+ FTP site](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).

### Installing Python dependencies

```bash
pip install -r requirements.txt
# optional, for .xlsx export:
pip install openpyxl
```

## Running the app

```bash
python main.py
```

## Using the app

1. **Query FASTA** — browse to your query protein/CDS FASTA (e.g. a DEG set).
2. **TF Database FASTA** — browse to your custom, family-annotated TF FASTA
   (PlantTFDB header format required — see above).
3. Adjust **BLAST Parameters** if needed:
   - **Program** — leave on *Auto-detect* to automatically pick `blastp` vs
     `blastx` based on the query sequence type.
   - **E-value threshold** — default `1e-5`.
   - **Max hits / query** — how many distinct database hits to retain per
     query; the single best hit determines the predicted TF family, and the
     rest are viewable by double-clicking a result row.
   - **Threads** — number of CPU threads BLAST should use.
   - **Force rebuild database** — check this if you edited the database
     FASTA in place without renaming it (the app normally caches built
     BLAST databases keyed by file path + size + modification time).
4. Click **Run BLAST Search**. Progress, live log messages, and a cancel
   button are shown while the search runs.
5. Once finished, use the **search box**, **family filter**, and
   **"Classified only"** checkbox above the results table to explore hits.
   Double-click any row to see all alternative hits retained for that query.
6. Use **Export Current View...** to save exactly what's currently
   filtered/sorted in the table, or **File → Save All Results As...** to
   save the complete, unfiltered result set. Both support `.tsv`, `.csv`,
   and `.xlsx` (if `openpyxl` is installed).

A **TF Family Distribution** tab below the results table shows a bar chart
summarizing how many query sequences were assigned to each TF family.

## Project layout

```
planttf_blast_gui/
├── main.py                 # entry point
├── requirements.txt
├── README.md
└── app/
    ├── __init__.py
    ├── blast_engine.py      # non-GUI: FASTA parsing, BLAST db caching, running BLAST
    ├── blast_worker.py      # QThread wrapper around blast_engine
    ├── results_model.py     # Qt table model + filter/sort proxy
    ├── family_chart.py      # dependency-free bar chart widget
    └── main_window.py       # the PyQt6 main window / UI
```

`app/blast_engine.py` has no PyQt dependency and can be reused or unit
tested standalone.
