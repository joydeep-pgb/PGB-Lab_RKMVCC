"""
blast_worker.py
----------------
QThread-based worker that runs the whole pipeline (build/reuse database,
run BLAST in chunks, parse + assemble results) off the GUI thread, emitting
Qt signals so MainWindow can update progress bars, logs and the results
table without blocking.
"""

from PyQt6.QtCore import QThread, pyqtSignal

from . import blast_engine as be


class BlastWorker(QThread):
    log_message = pyqtSignal(str)
    progress = pyqtSignal(int, int, str)  # done, total, message
    finished_ok = pyqtSignal(list, list)  # results (List[QueryResult]), family_dist
    failed = pyqtSignal(str)
    cancelled = pyqtSignal()

    def __init__(
        self,
        query_path: str,
        db_fasta_path: str,
        program: str,
        evalue: float,
        max_target_seqs: int,
        threads: int,
        force_rebuild_db: bool,
        chunk_size: int = 200,
        parent=None,
    ):
        super().__init__(parent)
        self.query_path = query_path
        self.db_fasta_path = db_fasta_path
        self.program = program
        self.evalue = evalue
        self.max_target_seqs = max_target_seqs
        self.threads = threads
        self.force_rebuild_db = force_rebuild_db
        self.chunk_size = chunk_size
        self._stop_requested = False

    def request_stop(self):
        self._stop_requested = True

    def _should_stop(self) -> bool:
        return self._stop_requested

    def run(self):
        try:
            if not be.blast_is_installed():
                self.failed.emit(be.install_hint())
                return

            self.log_message.emit("Parsing custom TF database headers...")
            db_entries = be.parse_all_db_headers(self.db_fasta_path)
            self.log_message.emit(
                f"Parsed {len(db_entries)} database sequences "
                f"({len(set(e.family for e in db_entries.values()))} distinct TF families)."
            )

            db_prefix = be.ensure_blast_db(
                self.db_fasta_path,
                log=lambda m: self.log_message.emit(m),
                force_rebuild=self.force_rebuild_db,
            )

            if self._should_stop():
                self.cancelled.emit()
                return

            self.log_message.emit(
                f"Running {self.program} (e-value <= {self.evalue:g}, "
                f"max hits/query = {self.max_target_seqs}, threads = {self.threads})..."
            )

            rows = be.run_blast_chunked(
                self.query_path,
                db_prefix,
                program=self.program,
                evalue=self.evalue,
                max_target_seqs=self.max_target_seqs,
                threads=self.threads,
                chunk_size=self.chunk_size,
                progress_cb=lambda d, t, m: self.progress.emit(d, t, m),
                should_stop=self._should_stop,
                log=lambda m: self.log_message.emit(m),
            )

            self.log_message.emit("Assembling PlantTFDB-style results...")
            results = be.build_query_results(
                self.query_path, rows, db_entries, max_hits_per_query=self.max_target_seqs
            )
            fam_dist = be.family_distribution(results)

            matched = sum(1 for r in results if r.matched)
            self.log_message.emit(
                f"Done. {matched}/{len(results)} query sequences classified into a TF family."
            )
            self.finished_ok.emit(results, fam_dist)

        except be.BlastCancelled:
            self.cancelled.emit()
        except be.BlastError as exc:
            self.failed.emit(str(exc))
        except Exception as exc:  # pragma: no cover - safety net
            self.failed.emit(f"Unexpected error: {exc}")
