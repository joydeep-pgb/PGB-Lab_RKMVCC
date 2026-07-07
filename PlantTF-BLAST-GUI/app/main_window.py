"""
main_window.py
----------------
The main PyQt6 window for the PlantTF-BLAST GUI application.
"""

import csv
import os
from typing import List, Optional

from PyQt6.QtCore import Qt
from PyQt6.QtGui import QAction, QColor
from PyQt6.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QComboBox,
    QDialog,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QPlainTextEdit,
    QProgressBar,
    QPushButton,
    QSizePolicy,
    QSpinBox,
    QSplitter,
    QStatusBar,
    QTableView,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from . import blast_engine as be
from .blast_worker import BlastWorker
from .family_chart import FamilyChartScrollArea
from .results_model import COLUMNS, ResultsFilterProxyModel, ResultsTableModel

APP_TITLE = "PlantTF-BLAST — Custom Transcription Factor BLAST Search"

EVALUE_PRESETS = ["1e-3", "1e-5", "1e-10", "1e-20", "1e-50", "1e-100"]


class HitDetailDialog(QDialog):
    """Shows all retained hits for a single query (alternative candidate families)."""

    def __init__(self, result, parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"Alternative hits — {result.query_id}")
        self.resize(720, 320)
        layout = QVBoxLayout(self)

        info = QLabel(
            f"<b>{result.query_id}</b> — {len(result.all_hits)} hit(s) retained "
            f"(best hit determines the predicted TF family)."
        )
        info.setWordWrap(True)
        layout.addWidget(info)

        table = QTableView(self)
        model = ResultsTableModel()
        # Wrap all_hits as pseudo QueryResult rows for display re-use
        from .blast_engine import QueryResult

        rows = []
        for h in result.all_hits:
            rows.append(
                QueryResult(
                    query_id=result.query_id,
                    matched=True,
                    family=h.family,
                    hit_id=h.hit_id,
                    evalue=h.evalue,
                    bitscore=h.bitscore,
                    pident=h.pident,
                    qcov=h.qcov,
                    description=h.description,
                    species=h.species,
                )
            )
        model.set_results(rows)
        table.setModel(model)
        table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeMode.ResizeToContents
        )
        table.horizontalHeader().setStretchLastSection(True)
        table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        layout.addWidget(table)

        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.accept)
        layout.addWidget(close_btn, alignment=Qt.AlignmentFlag.AlignRight)


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle(APP_TITLE)
        self.resize(1180, 720)

        self.query_path: Optional[str] = None
        self.db_path: Optional[str] = None
        self.results: List = []
        self.family_dist = []
        self.worker: Optional[BlastWorker] = None

        self._build_ui()
        self._check_blast_on_startup()

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------
    def _build_ui(self):
        self._build_menu()

        central = QWidget()
        self.setCentralWidget(central)
        root_layout = QHBoxLayout(central)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        root_layout.addWidget(splitter)

        splitter.addWidget(self._build_left_panel())
        splitter.addWidget(self._build_right_panel())
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([360, 820])

        self.setStatusBar(QStatusBar())
        self.statusBar().showMessage("Ready.")

    def _build_menu(self):
        menu = self.menuBar()

        file_menu = menu.addMenu("&File")
        act_open_query = QAction("Open &Query FASTA...", self)
        act_open_query.triggered.connect(self.browse_query)
        act_open_db = QAction("Open TF &Database FASTA...", self)
        act_open_db.triggered.connect(self.browse_db)
        act_save = QAction("&Save All Results As...", self)
        act_save.triggered.connect(self.export_all_results)
        act_exit = QAction("E&xit", self)
        act_exit.triggered.connect(self.close)
        file_menu.addAction(act_open_query)
        file_menu.addAction(act_open_db)
        file_menu.addSeparator()
        file_menu.addAction(act_save)
        file_menu.addSeparator()
        file_menu.addAction(act_exit)

        run_menu = menu.addMenu("&Run")
        self.act_run = QAction("Run BLAST Search", self)
        self.act_run.triggered.connect(self.run_blast)
        self.act_cancel = QAction("Cancel", self)
        self.act_cancel.triggered.connect(self.cancel_blast)
        self.act_cancel.setEnabled(False)
        run_menu.addAction(self.act_run)
        run_menu.addAction(self.act_cancel)

        help_menu = menu.addMenu("&Help")
        act_how = QAction("How This Works", self)
        act_how.triggered.connect(self.show_how_it_works)
        act_about = QAction("About", self)
        act_about.triggered.connect(self.show_about)
        help_menu.addAction(act_how)
        help_menu.addAction(act_about)

    def _build_left_panel(self) -> QWidget:
        panel = QWidget()
        layout = QVBoxLayout(panel)
        layout.setAlignment(Qt.AlignmentFlag.AlignTop)

        # --- Input files ---
        files_box = QGroupBox("Input Files")
        files_layout = QFormLayout()
        files_box.setLayout(files_layout)

        self.query_edit = QLineEdit()
        self.query_edit.setReadOnly(True)
        self.query_edit.setPlaceholderText("No query FASTA selected")
        query_btn = QPushButton("Browse...")
        query_btn.clicked.connect(self.browse_query)
        query_row = QHBoxLayout()
        query_row.addWidget(self.query_edit)
        query_row.addWidget(query_btn)
        files_layout.addRow("Query FASTA:", self._wrap(query_row))

        self.db_edit = QLineEdit()
        self.db_edit.setReadOnly(True)
        self.db_edit.setPlaceholderText("No custom TF database selected")
        db_btn = QPushButton("Browse...")
        db_btn.clicked.connect(self.browse_db)
        db_row = QHBoxLayout()
        db_row.addWidget(self.db_edit)
        db_row.addWidget(db_btn)
        files_layout.addRow("TF Database FASTA:", self._wrap(db_row))

        self.files_info_label = QLabel(
            "Database headers must follow the PlantTFDB convention:\n"
            ">GeneID Species|TF_Family|Description"
        )
        self.files_info_label.setWordWrap(True)
        self.files_info_label.setStyleSheet("color: #666666; font-size: 11px;")
        files_layout.addRow(self.files_info_label)

        layout.addWidget(files_box)

        # --- Parameters ---
        params_box = QGroupBox("BLAST Parameters")
        params_layout = QFormLayout()
        params_box.setLayout(params_layout)

        self.program_combo = QComboBox()
        self.program_combo.addItems(["Auto-detect", "blastp (protein query)", "blastx (nucleotide query)"])
        params_layout.addRow("Program:", self.program_combo)

        self.evalue_combo = QComboBox()
        self.evalue_combo.setEditable(True)
        self.evalue_combo.addItems(EVALUE_PRESETS)
        self.evalue_combo.setCurrentText("1e-5")
        params_layout.addRow("E-value threshold:", self.evalue_combo)

        self.max_hits_spin = QSpinBox()
        self.max_hits_spin.setRange(1, 20)
        self.max_hits_spin.setValue(5)
        self.max_hits_spin.setToolTip(
            "Number of distinct database hits retained per query.\n"
            "The single best hit determines the predicted TF family."
        )
        params_layout.addRow("Max hits / query:", self.max_hits_spin)

        self.threads_spin = QSpinBox()
        self.threads_spin.setRange(1, max(1, os.cpu_count() or 1))
        self.threads_spin.setValue(min(4, os.cpu_count() or 1))
        params_layout.addRow("Threads:", self.threads_spin)

        self.force_rebuild_check = QCheckBox("Force rebuild database")
        self.force_rebuild_check.setToolTip(
            "Rebuild the BLAST database even if a cached copy exists\n"
            "(use this if the database FASTA content changed but the filename/path did not)."
        )
        params_layout.addRow(self.force_rebuild_check)

        layout.addWidget(params_box)

        # --- Run controls ---
        run_box = QGroupBox("Run")
        run_layout = QVBoxLayout()
        run_box.setLayout(run_layout)

        btn_row = QHBoxLayout()
        self.run_btn = QPushButton("Run BLAST Search")
        self.run_btn.setStyleSheet("font-weight: bold; padding: 6px;")
        self.run_btn.clicked.connect(self.run_blast)
        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.clicked.connect(self.cancel_blast)
        self.cancel_btn.setEnabled(False)
        btn_row.addWidget(self.run_btn)
        btn_row.addWidget(self.cancel_btn)
        run_layout.addLayout(btn_row)

        self.progress_bar = QProgressBar()
        self.progress_bar.setValue(0)
        run_layout.addWidget(self.progress_bar)

        self.run_status_label = QLabel("")
        self.run_status_label.setWordWrap(True)
        self.run_status_label.setStyleSheet("color: #666666; font-size: 11px;")
        run_layout.addWidget(self.run_status_label)

        layout.addWidget(run_box)

        # --- Summary ---
        summary_box = QGroupBox("Summary")
        summary_layout = QFormLayout()
        summary_box.setLayout(summary_layout)
        self.summary_total = QLabel("-")
        self.summary_matched = QLabel("-")
        self.summary_unmatched = QLabel("-")
        self.summary_families = QLabel("-")
        summary_layout.addRow("Total queries:", self.summary_total)
        summary_layout.addRow("Classified as TF:", self.summary_matched)
        summary_layout.addRow("No hit:", self.summary_unmatched)
        summary_layout.addRow("Distinct families found:", self.summary_families)
        layout.addWidget(summary_box)

        layout.addStretch(1)
        return panel

    def _build_right_panel(self) -> QWidget:
        panel = QWidget()
        layout = QVBoxLayout(panel)

        # toolbar row: search / filter / export
        toolbar = QHBoxLayout()
        self.search_edit = QLineEdit()
        self.search_edit.setPlaceholderText("Search query ID, family, hit ID, description...")
        self.search_edit.textChanged.connect(self._apply_filters)
        toolbar.addWidget(self.search_edit, stretch=2)

        self.family_filter_combo = QComboBox()
        self.family_filter_combo.addItem("All families")
        self.family_filter_combo.currentTextChanged.connect(self._apply_filters)
        toolbar.addWidget(self.family_filter_combo, stretch=1)

        self.matched_only_check = QCheckBox("Classified only")
        self.matched_only_check.stateChanged.connect(self._apply_filters)
        toolbar.addWidget(self.matched_only_check)

        export_btn = QPushButton("Export Current View...")
        export_btn.clicked.connect(self.export_filtered_results)
        toolbar.addWidget(export_btn)

        layout.addLayout(toolbar)

        # results table
        self.model = ResultsTableModel()
        self.proxy = ResultsFilterProxyModel()
        self.proxy.setSourceModel(self.model)

        self.table = QTableView()
        self.table.setModel(self.proxy)
        self.table.setSortingEnabled(True)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeMode.Interactive
        )
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setAlternatingRowColors(True)
        self.table.doubleClicked.connect(self._show_hit_detail)
        for i, (name, _) in enumerate(COLUMNS):
            width = 220 if name == "Description" else 110
            self.table.setColumnWidth(i, width)

        bottom_tabs = QTabWidget()
        self.family_chart = FamilyChartScrollArea()
        bottom_tabs.addTab(self.family_chart, "TF Family Distribution")

        self.log_view = QPlainTextEdit()
        self.log_view.setReadOnly(True)
        self.log_view.setMaximumBlockCount(5000)
        bottom_tabs.addTab(self.log_view, "Log")

        vsplit = QSplitter(Qt.Orientation.Vertical)
        vsplit.addWidget(self.table)
        vsplit.addWidget(bottom_tabs)
        vsplit.setStretchFactor(0, 3)
        vsplit.setStretchFactor(1, 2)

        layout.addWidget(vsplit)
        return panel

    @staticmethod
    def _wrap(inner_layout) -> QWidget:
        w = QWidget()
        w.setLayout(inner_layout)
        return w

    # ------------------------------------------------------------------
    # Startup checks
    # ------------------------------------------------------------------
    def _check_blast_on_startup(self):
        if not be.blast_is_installed():
            QMessageBox.warning(self, "NCBI BLAST+ not found", be.install_hint())
            self.statusBar().showMessage("NCBI BLAST+ not found on PATH.")

    # ------------------------------------------------------------------
    # File selection
    # ------------------------------------------------------------------
    def browse_query(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select Query FASTA File", "", "FASTA files (*.fasta *.fa *.faa *.pep *.txt);;All files (*)"
        )
        if not path:
            return
        self.query_path = path
        self.query_edit.setText(path)
        try:
            n = be.count_fasta_sequences(path)
            seq_type = be.detect_sequence_type(path)
            self.statusBar().showMessage(
                f"Query loaded: {os.path.basename(path)} — {n} sequences, detected as {seq_type}."
            )
        except Exception as exc:
            QMessageBox.critical(self, "Error reading query file", str(exc))

    def browse_db(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select Custom TF Database FASTA File", "", "FASTA files (*.fasta *.fa *.faa *.pep *.txt);;All files (*)"
        )
        if not path:
            return
        self.db_path = path
        self.db_edit.setText(path)
        try:
            entries = be.parse_all_db_headers(path)
            n_fam = len(set(e.family for e in entries.values()))
            self.statusBar().showMessage(
                f"Database loaded: {os.path.basename(path)} — {len(entries)} sequences, "
                f"{n_fam} distinct TF families detected from headers."
            )
            families = sorted(set(e.family for e in entries.values()))
            self.family_filter_combo.blockSignals(True)
            self.family_filter_combo.clear()
            self.family_filter_combo.addItem("All families")
            self.family_filter_combo.addItems(families)
            self.family_filter_combo.blockSignals(False)
        except Exception as exc:
            QMessageBox.critical(self, "Error reading database file", str(exc))

    # ------------------------------------------------------------------
    # Running BLAST
    # ------------------------------------------------------------------
    def _resolve_program(self) -> str:
        choice = self.program_combo.currentText()
        if choice.startswith("blastp"):
            return "blastp"
        if choice.startswith("blastx"):
            return "blastx"
        # auto-detect
        seq_type = be.detect_sequence_type(self.query_path)
        return "blastp" if seq_type == "protein" else "blastx"

    def _parse_evalue(self) -> Optional[float]:
        text = self.evalue_combo.currentText().strip()
        try:
            return float(text)
        except ValueError:
            return None

    def run_blast(self):
        if self.worker is not None and self.worker.isRunning():
            QMessageBox.information(self, "Already running", "A BLAST search is already running.")
            return

        if not self.query_path or not os.path.isfile(self.query_path):
            QMessageBox.warning(self, "Missing query", "Please select a valid query FASTA file first.")
            return
        if not self.db_path or not os.path.isfile(self.db_path):
            QMessageBox.warning(self, "Missing database", "Please select a valid custom TF database FASTA file first.")
            return
        if not be.blast_is_installed():
            QMessageBox.critical(self, "NCBI BLAST+ not found", be.install_hint())
            return

        evalue = self._parse_evalue()
        if evalue is None or evalue <= 0:
            QMessageBox.warning(
                self, "Invalid e-value", "Please enter a valid e-value threshold, e.g. 1e-5."
            )
            return

        program = self._resolve_program()

        self.log_view.clear()
        self.progress_bar.setValue(0)
        self.run_status_label.setText("Starting...")
        self._set_running_state(True)

        self.worker = BlastWorker(
            query_path=self.query_path,
            db_fasta_path=self.db_path,
            program=program,
            evalue=evalue,
            max_target_seqs=self.max_hits_spin.value(),
            threads=self.threads_spin.value(),
            force_rebuild_db=self.force_rebuild_check.isChecked(),
        )
        self.worker.log_message.connect(self._on_log)
        self.worker.progress.connect(self._on_progress)
        self.worker.finished_ok.connect(self._on_finished)
        self.worker.failed.connect(self._on_failed)
        self.worker.cancelled.connect(self._on_cancelled)
        self.worker.start()

    def cancel_blast(self):
        if self.worker is not None and self.worker.isRunning():
            self.worker.request_stop()
            self.run_status_label.setText("Cancelling...")

    def _set_running_state(self, running: bool):
        self.run_btn.setEnabled(not running)
        self.act_run.setEnabled(not running)
        self.cancel_btn.setEnabled(running)
        self.act_cancel.setEnabled(running)
        for w in (
            self.query_edit,
            self.db_edit,
            self.program_combo,
            self.evalue_combo,
            self.max_hits_spin,
            self.threads_spin,
            self.force_rebuild_check,
        ):
            w.setEnabled(not running)

    def _on_log(self, message: str):
        self.log_view.appendPlainText(message)

    def _on_progress(self, done: int, total: int, message: str):
        pct = int(100 * done / max(1, total))
        self.progress_bar.setValue(pct)
        self.run_status_label.setText(message)
        self.statusBar().showMessage(message)

    def _on_finished(self, results, family_dist):
        self.results = results
        self.family_dist = family_dist
        self.model.set_results(results)

        families = sorted(set(r.family for r in results if r.matched))
        current = self.family_filter_combo.currentText()
        self.family_filter_combo.blockSignals(True)
        self.family_filter_combo.clear()
        self.family_filter_combo.addItem("All families")
        self.family_filter_combo.addItems(families)
        idx = self.family_filter_combo.findText(current)
        self.family_filter_combo.setCurrentIndex(idx if idx >= 0 else 0)
        self.family_filter_combo.blockSignals(False)

        self.family_chart.set_data(family_dist)

        total = len(results)
        matched = sum(1 for r in results if r.matched)
        self.summary_total.setText(str(total))
        self.summary_matched.setText(str(matched))
        self.summary_unmatched.setText(str(total - matched))
        self.summary_families.setText(str(len(families)))

        self.progress_bar.setValue(100)
        self.run_status_label.setText("Finished.")
        self.statusBar().showMessage(
            f"Finished: {matched}/{total} query sequences classified into a TF family."
        )
        self._set_running_state(False)

    def _on_failed(self, message: str):
        self._set_running_state(False)
        self.run_status_label.setText("Failed.")
        QMessageBox.critical(self, "BLAST run failed", message)

    def _on_cancelled(self):
        self._set_running_state(False)
        self.run_status_label.setText("Cancelled.")
        self.statusBar().showMessage("BLAST run cancelled.")

    # ------------------------------------------------------------------
    # Filtering
    # ------------------------------------------------------------------
    def _apply_filters(self):
        self.proxy.set_search_text(self.search_edit.text())
        self.proxy.set_family_filter(self.family_filter_combo.currentText())
        self.proxy.set_matched_only(self.matched_only_check.isChecked())

    def _show_hit_detail(self, index):
        source_index = self.proxy.mapToSource(index)
        result = self.model.result_at(source_index.row())
        if result is None or not result.matched or len(result.all_hits) <= 1:
            return
        dlg = HitDetailDialog(result, self)
        dlg.exec()

    # ------------------------------------------------------------------
    # Export
    # ------------------------------------------------------------------
    def _write_results(self, path: str, results: List):
        headers = [name for name, _ in COLUMNS]
        is_excel = path.lower().endswith((".xlsx", ".xls"))
        rows = []
        for r in results:
            if not r.matched:
                rows.append([r.query_id, "No hit", "-", "-", "-", "-", "-"])
                continue
            evalue_str = f"{r.evalue:.2e}" if r.evalue else "-"
            pident_str = f"{r.pident:.1f}" if r.pident is not None else "-"
            qcov_str = f"{r.qcov:.1f}" if r.qcov is not None else "-"
            rows.append(
                [r.query_id, r.family, r.hit_id, evalue_str, pident_str, qcov_str, r.description]
            )

        if is_excel:
            try:
                from openpyxl import Workbook
            except ImportError:
                raise RuntimeError(
                    "Saving as .xlsx requires the 'openpyxl' package "
                    "(pip install openpyxl). Please choose .tsv or .csv instead."
                )
            wb = Workbook()
            ws = wb.active
            ws.title = "TF_BLAST_results"
            ws.append(headers)
            for row in rows:
                ws.append(row)
            wb.save(path)
        else:
            delimiter = "," if path.lower().endswith(".csv") else "\t"
            with open(path, "w", newline="") as fh:
                writer = csv.writer(fh, delimiter=delimiter)
                writer.writerow(headers)
                writer.writerows(rows)

    def export_all_results(self):
        if not self.results:
            QMessageBox.information(self, "Nothing to export", "Run a BLAST search first.")
            return
        self._export_dialog(self.results)

    def export_filtered_results(self):
        if not self.results:
            QMessageBox.information(self, "Nothing to export", "Run a BLAST search first.")
            return
        visible_results = []
        for row in range(self.proxy.rowCount()):
            source_index = self.proxy.mapToSource(self.proxy.index(row, 0))
            r = self.model.result_at(source_index.row())
            if r is not None:
                visible_results.append(r)
        self._export_dialog(visible_results)

    def _export_dialog(self, results: List):
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Save Results",
            "tf_blast_results.tsv",
            "Tab-separated (*.tsv);;CSV (*.csv);;Excel (*.xlsx)",
        )
        if not path:
            return
        try:
            self._write_results(path, results)
            self.statusBar().showMessage(f"Saved {len(results)} rows to {path}")
        except Exception as exc:
            QMessageBox.critical(self, "Export failed", str(exc))

    # ------------------------------------------------------------------
    # Help
    # ------------------------------------------------------------------
    def show_how_it_works(self):
        QMessageBox.information(
            self,
            "How This Works",
            "<b>PlantTF-BLAST</b> predicts transcription factor (TF) identity for your "
            "query protein sequences by homology search against a custom, "
            "family-annotated TF database (e.g. a PlantTFDB species download).<br><br>"
            "1. Your TF database FASTA headers must follow the format:<br>"
            "&nbsp;&nbsp;<code>&gt;GeneID Species|TF_Family|Description</code><br><br>"
            "2. A local BLAST protein database is built (and cached) from that FASTA.<br><br>"
            "3. Your query sequences are searched (blastp, or blastx for nucleotide "
            "queries) against that database.<br><br>"
            "4. For each query, the best-scoring hit's TF family is reported as the "
            "<i>predicted TF family</i>, along with the hit ID, e-value and description "
            "— the same style of output produced by PlantTFDB's BLAST tool.<br><br>"
            "Note: this is a homology-based (BLAST) classification, not a Pfam/HMM "
            "domain-based classification like PlantTFDB's own internal pipeline; "
            "results depend entirely on the quality and coverage of your custom "
            "database's family annotations.",
        )

    def show_about(self):
        QMessageBox.about(
            self,
            "About PlantTF-BLAST",
            "<b>PlantTF-BLAST GUI</b><br>"
            "A local, custom-database BLAST tool for predicting transcription factor "
            "families in plant proteomes, styled after PlantTFDB's BLAST search "
            "results.<br><br>"
            "Built with PyQt6 and NCBI BLAST+.",
        )

    # ------------------------------------------------------------------
    def closeEvent(self, event):
        if self.worker is not None and self.worker.isRunning():
            self.worker.request_stop()
            self.worker.wait(3000)
        event.accept()
