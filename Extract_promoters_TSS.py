#!/usr/bin/env python3
"""
Promoter Sequence Extractor
----------------------------
A PyQt5 GUI for extracting promoter sequences of selected genes from a
Phytozome-style GFF3 annotation, using the longest transcript per gene.

Dependencies:
    pip install PyQt5 biopython
    (or: conda install -c conda-forge -c bioconda pyqt biopython)

Run:
    python promoter_extractor_gui.py
"""

import os
import re
import sys
import traceback

from PyQt5.QtCore import Qt, QThread, pyqtSignal, QUrl
from PyQt5.QtGui import QFont, QDesktopServices, QTextCursor
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QFormLayout,
    QGroupBox, QLabel, QLineEdit, QPushButton, QPlainTextEdit, QTextEdit,
    QSpinBox, QFileDialog, QProgressBar, QMessageBox, QSplitter, QStatusBar,
    QAction, QSizePolicy, QCheckBox
)

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

APP_TITLE = "Promoter Sequence Extractor"
APP_VERSION = "1.0"


# --------------------------------------------------------------------------
# Worker thread — keeps the extraction logic off the GUI thread
# --------------------------------------------------------------------------
class ExtractionWorker(QThread):
    log = pyqtSignal(str)
    progress = pyqtSignal(int, int)          # current, total  (-1, -1 = busy/indeterminate)
    finished_ok = pyqtSignal(int, list, str)  # extracted_count, missing_genes, output_path
    failed = pyqtSignal(str)

    VERSION_SUFFIX_RE = re.compile(r"\.v\d+(?:\.\d+)?$")

    def __init__(self, genome_path, gff_path, gene_ids, output_path, promoter_length,
                 strip_version=True):
        super().__init__()
        self.genome_path = genome_path
        self.gff_path = gff_path
        self.gene_ids = gene_ids
        self.output_path = output_path
        self.promoter_length = promoter_length
        self.strip_version = strip_version
        self._cancelled = False

    def cancel(self):
        self._cancelled = True

    def run(self):
        try:
            self.progress.emit(-1, -1)
            self.log.emit("Loading genome FASTA ...")
            genome = SeqIO.to_dict(SeqIO.parse(self.genome_path, "fasta"))
            self.log.emit(f"Loaded {len(genome)} sequence(s) from genome.")

            target_genes = set(self.gene_ids)
            self.log.emit(f"Target genes: {len(target_genes)}")

            records = []
            missing_chroms = set()

            # Diagnostics — help pinpoint *why* nothing matched, rather than
            # just reporting every gene as "not found".
            total_mRNA = 0
            with_longest_tag = 0
            parent_id_samples = []      # Parent= values seen on longest=1 mRNAs
            attr_format_sample = None   # one raw attribute string, for reference

            total_size = max(1, os.path.getsize(self.gff_path))
            self.log.emit("Scanning annotation (GFF3) ...")

            with open(self.gff_path) as gff:
                i = 0
                while True:
                    line = gff.readline()
                    if not line:
                        break
                    i += 1

                    if self._cancelled:
                        self.log.emit("Cancelled by user.")
                        return

                    if i % 2000 == 0:
                        self.progress.emit(gff.tell(), total_size)

                    if line.startswith("#"):
                        continue
                    cols = line.rstrip("\n").split("\t")
                    if len(cols) != 9:
                        continue
                    chrom, source, feature, start, end, score, strand, phase, attr = cols
                    if feature != "mRNA":
                        continue
                    total_mRNA += 1
                    if attr_format_sample is None:
                        attr_format_sample = attr
                    if "longest=1" not in attr:
                        continue
                    with_longest_tag += 1

                    attrs = {}
                    for item in attr.split(";"):
                        if "=" in item:
                            k, v = item.split("=", 1)
                            attrs[k] = v

                    raw_parent = attrs.get("Parent")
                    if len(parent_id_samples) < 5 and raw_parent not in parent_id_samples:
                        parent_id_samples.append(raw_parent)

                    if self.strip_version and raw_parent:
                        gene_id = self.VERSION_SUFFIX_RE.sub("", raw_parent)
                    else:
                        gene_id = raw_parent

                    if gene_id not in target_genes:
                        continue

                    transcript_id = attrs.get("ID")
                    start = int(start)
                    end = int(end)

                    if chrom not in genome:
                        missing_chroms.add(chrom)
                        self.log.emit(
                            f"WARNING: '{chrom}' (gene {gene_id}) not found in genome FASTA — skipped."
                        )
                        continue

                    seq = genome[chrom].seq
                    chr_len = len(seq)

                    if strand == "+":
                        tss = start
                        p_start = max(0, tss - self.promoter_length - 1)
                        p_end = tss - 1
                        promoter = seq[p_start:p_end]
                    else:
                        tss = end
                        p_start = tss
                        p_end = min(chr_len, tss + self.promoter_length)
                        promoter = seq[p_start:p_end].reverse_complement()

                    header = (
                        f"{gene_id}"
                        f"|{transcript_id}"
                        f"|{chrom}"
                        f"|{strand}"
                        f"|TSS={tss}"
                        f"|Promoter={len(promoter)}bp"
                    )
                    records.append(SeqRecord(promoter, id=header, description=""))
                    self.log.emit(f"Extracted {gene_id} ({transcript_id}) — {len(promoter)} bp")

            self.progress.emit(total_size, total_size)

            SeqIO.write(records, self.output_path, "fasta")

            found = {r.id.split("|")[0] for r in records}
            missing = sorted(target_genes - found)

            # --- Diagnostics summary -------------------------------------
            self.log.emit("\n--- Diagnostics ---")
            self.log.emit(f"mRNA features scanned: {total_mRNA}")
            self.log.emit(f"...with a 'longest=1' tag: {with_longest_tag}")
            self.log.emit(f"Matched to your gene ID list: {len(records)}")

            if total_mRNA == 0:
                self.log.emit(
                    "No 'mRNA' features were found at all in this GFF3. "
                    "Check that column 3 (feature type) uses exactly 'mRNA' "
                    "(some GFF3s use 'transcript' instead)."
                )
            elif with_longest_tag == 0:
                self.log.emit(
                    "No mRNA line contains the literal tag 'longest=1'. "
                    "This GFF3 may mark the primary/representative transcript "
                    "differently. Example attribute string from this file:\n"
                    f"  {attr_format_sample}"
                )
            elif len(records) == 0 and parent_id_samples:
                sample_targets = list(target_genes)[:5]
                if self.strip_version:
                    stripped_samples = [
                        self.VERSION_SUFFIX_RE.sub("", p) if p else p
                        for p in parent_id_samples
                    ]
                    self.log.emit(
                        "Found candidate transcripts, but none matched your gene "
                        "ID list even after stripping '.vX.Y' version suffixes. "
                        "Compare:\n"
                        f"  Parent= IDs in GFF3 (raw sample):     {parent_id_samples}\n"
                        f"  Parent= IDs after stripping (sample): {stripped_samples}\n"
                        f"  Your gene IDs (sample):               {sample_targets}"
                    )
                else:
                    self.log.emit(
                        "Found candidate transcripts, but none matched your gene "
                        "ID list — likely an ID-formatting mismatch (version "
                        "suffix, case, extra whitespace, etc.). Try enabling "
                        "'Strip Phytozome version suffix' above. Compare:\n"
                        f"  Parent= IDs in GFF3 (sample): {parent_id_samples}\n"
                        f"  Your gene IDs (sample):       {sample_targets}"
                    )

            self.finished_ok.emit(len(records), missing, self.output_path)

        except Exception:
            self.failed.emit(traceback.format_exc())


# --------------------------------------------------------------------------
# Main window
# --------------------------------------------------------------------------
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.worker = None
        self.setWindowTitle(f"{APP_TITLE}")
        self.resize(920, 760)
        self._build_ui()
        self._apply_style()

    # ---- UI construction -------------------------------------------------
    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(14, 14, 14, 14)
        root.setSpacing(12)

        # --- Menu bar ---
        file_menu = self.menuBar().addMenu("&File")
        exit_action = QAction("Exit", self)
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)

        help_menu = self.menuBar().addMenu("&Help")
        about_action = QAction("About", self)
        about_action.triggered.connect(self._show_about)
        help_menu.addAction(about_action)

        # --- Inputs group ---
        inputs_group = QGroupBox("Input Files")
        form = QFormLayout()
        form.setSpacing(10)

        self.genome_edit = QLineEdit()
        self.genome_edit.setPlaceholderText("Path to genome FASTA file ...")
        genome_browse = QPushButton("Browse...")
        genome_browse.clicked.connect(self._browse_genome)
        form.addRow("Genome FASTA:", self._with_button(self.genome_edit, genome_browse))

        self.gff_edit = QLineEdit()
        self.gff_edit.setPlaceholderText("Path to Phytozome GFF3 annotation file ...")
        gff_browse = QPushButton("Browse...")
        gff_browse.clicked.connect(self._browse_gff)
        form.addRow("GFF3 Annotation:", self._with_button(self.gff_edit, gff_browse))

        self.length_spin = QSpinBox()
        self.length_spin.setRange(1, 100000)
        self.length_spin.setValue(1500)
        self.length_spin.setSuffix(" bp")
        self.length_spin.setMaximumWidth(140)
        form.addRow("Promoter length:", self.length_spin)

        self.strip_version_check = QCheckBox(
            "Strip Phytozome version suffix (e.g. \".v2.1\") from gene IDs when matching"
        )
        self.strip_version_check.setChecked(True)
        self.strip_version_check.setToolTip(
            "Phytozome GFF3 files often label transcripts as Parent=GeneID.v2.1.\n"
            "Enable this if your gene ID list uses the bare ID (GeneID) without the version."
        )
        form.addRow("", self.strip_version_check)

        inputs_group.setLayout(form)
        root.addWidget(inputs_group)

        # --- Gene ID list group ---
        genes_group = QGroupBox("Gene ID List")
        genes_layout = QVBoxLayout()

        genes_header = QHBoxLayout()
        genes_header.addWidget(QLabel("Paste gene IDs below (one per line, or separated by commas/spaces):"))
        genes_header.addStretch()
        load_ids_btn = QPushButton("Load from file...")
        load_ids_btn.clicked.connect(self._load_ids_from_file)
        clear_ids_btn = QPushButton("Clear")
        clear_ids_btn.clicked.connect(lambda: self.ids_text.clear())
        genes_header.addWidget(load_ids_btn)
        genes_header.addWidget(clear_ids_btn)
        genes_layout.addLayout(genes_header)

        self.ids_text = QPlainTextEdit()
        self.ids_text.setPlaceholderText("Phvul.008G220400\nPhvul.001G045300\n...")
        self.ids_text.setFont(QFont("Monospace", 10))
        self.ids_text.setFixedHeight(130)
        genes_layout.addWidget(self.ids_text)

        self.gene_count_label = QLabel("0 gene ID(s) entered")
        self.gene_count_label.setObjectName("hint")
        genes_layout.addWidget(self.gene_count_label)
        self.ids_text.textChanged.connect(self._update_gene_count)

        genes_group.setLayout(genes_layout)
        root.addWidget(genes_group)

        # --- Output group ---
        output_group = QGroupBox("Output")
        output_form = QFormLayout()
        self.output_edit = QLineEdit()
        self.output_edit.setPlaceholderText("Path to output FASTA file ...")
        output_browse = QPushButton("Save As...")
        output_browse.clicked.connect(self._browse_output)
        output_form.addRow("Output FASTA:", self._with_button(self.output_edit, output_browse))
        output_group.setLayout(output_form)
        root.addWidget(output_group)

        # --- Action buttons ---
        action_row = QHBoxLayout()
        self.run_btn = QPushButton("Extract Promoters")
        self.run_btn.setObjectName("primaryButton")
        self.run_btn.clicked.connect(self._start_extraction)

        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.setEnabled(False)
        self.cancel_btn.clicked.connect(self._cancel_extraction)

        self.open_folder_btn = QPushButton("Open Output Folder")
        self.open_folder_btn.setEnabled(False)
        self.open_folder_btn.clicked.connect(self._open_output_folder)

        action_row.addWidget(self.run_btn)
        action_row.addWidget(self.cancel_btn)
        action_row.addStretch()
        action_row.addWidget(self.open_folder_btn)
        root.addLayout(action_row)

        # --- Progress bar ---
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        self.progress_bar.setTextVisible(True)
        root.addWidget(self.progress_bar)

        # --- Log console ---
        log_group = QGroupBox("Log")
        log_layout = QVBoxLayout()
        self.log_view = QTextEdit()
        self.log_view.setReadOnly(True)
        self.log_view.setFont(QFont("Monospace", 9))
        log_layout.addWidget(self.log_view)
        log_group.setLayout(log_layout)
        root.addWidget(log_group, stretch=1)

        # --- Status bar ---
        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)
        self.status_bar.showMessage("Ready.")

    def _with_button(self, line_edit, button):
        w = QWidget()
        h = QHBoxLayout(w)
        h.setContentsMargins(0, 0, 0, 0)
        h.addWidget(line_edit)
        h.addWidget(button)
        return w

    def _apply_style(self):
        self.setStyleSheet("""
            QMainWindow { background-color: #f4f6f8; }
            QGroupBox {
                font-weight: 600;
                border: 1px solid #d0d7de;
                border-radius: 6px;
                margin-top: 10px;
                padding: 10px;
                background-color: #ffffff;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 4px;
                color: #2c3e50;
            }
            QLineEdit, QPlainTextEdit, QTextEdit, QSpinBox {
                border: 1px solid #c7cdd3;
                border-radius: 4px;
                padding: 4px 6px;
                background-color: #ffffff;
            }
            QPushButton {
                border: 1px solid #c7cdd3;
                border-radius: 4px;
                padding: 6px 14px;
                background-color: #ffffff;
            }
            QPushButton:hover { background-color: #eef2f6; }
            QPushButton:disabled { color: #9aa5b1; }
            QPushButton#primaryButton {
                background-color: #2f6feb;
                color: white;
                font-weight: 600;
                border: none;
            }
            QPushButton#primaryButton:hover { background-color: #2758c4; }
            QPushButton#primaryButton:disabled { background-color: #9db5ee; }
            QLabel#hint { color: #6a7480; font-size: 11px; }
            QProgressBar {
                border: 1px solid #c7cdd3;
                border-radius: 4px;
                text-align: center;
                background-color: #ffffff;
            }
            QProgressBar::chunk { background-color: #2f6feb; }
        """)

    # ---- Helpers -----------------------------------------------------
    def _parse_gene_ids(self):
        raw = self.ids_text.toPlainText()
        ids = set()
        for chunk in raw.replace(",", "\n").split("\n"):
            gid = chunk.strip()
            if gid:
                ids.add(gid)
        return ids

    def _update_gene_count(self):
        n = len(self._parse_gene_ids())
        self.gene_count_label.setText(f"{n} gene ID(s) entered")

    def _browse_genome(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select genome FASTA", "",
            "FASTA files (*.fa *.fasta *.fna *.fa.gz *.fasta.gz);;All files (*)"
        )
        if path:
            self.genome_edit.setText(path)

    def _browse_gff(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select GFF3 annotation", "",
            "GFF3 files (*.gff3 *.gff *.gff3.gz *.gff.gz);;All files (*)"
        )
        if path:
            self.gff_edit.setText(path)

    def _browse_output(self):
        path, _ = QFileDialog.getSaveFileName(
            self, "Save output FASTA as", "promoters.fasta",
            "FASTA files (*.fasta *.fa);;All files (*)"
        )
        if path:
            self.output_edit.setText(path)

    def _load_ids_from_file(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Load gene IDs from file", "", "Text files (*.txt);;All files (*)"
        )
        if not path:
            return
        try:
            with open(path) as f:
                content = f.read()
            current = self.ids_text.toPlainText()
            sep = "\n" if current.strip() else ""
            self.ids_text.setPlainText(current + sep + content.strip())
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Could not read file:\n{e}")

    def _open_output_folder(self):
        out_path = self.output_edit.text().strip()
        if out_path:
            folder = os.path.dirname(os.path.abspath(out_path)) or "."
            QDesktopServices.openUrl(QUrl.fromLocalFile(folder))

    def _show_about(self):
        QMessageBox.information(
            self, "About",
            f"<b>{APP_TITLE}</b> v{APP_VERSION}<br><br>"
            "Extracts upstream promoter sequences for a set of genes from a "
            "genome FASTA and a Phytozome-style GFF3 annotation, using the "
            "longest transcript per gene (attribute <i>longest=1</i>).<br><br>"
            "Built with PyQt5 and Biopython."
        )

    def _append_log(self, text):
        self.log_view.moveCursor(QTextCursor.End)
        self.log_view.insertPlainText(text + "\n")
        self.log_view.moveCursor(QTextCursor.End)

    def _set_inputs_enabled(self, enabled):
        for w in (self.genome_edit, self.gff_edit, self.length_spin,
                  self.ids_text, self.output_edit, self.run_btn):
            w.setEnabled(enabled)
        self.cancel_btn.setEnabled(not enabled)

    # ---- Validation & run ---------------------------------------------
    def _start_extraction(self):
        genome_path = self.genome_edit.text().strip()
        gff_path = self.gff_edit.text().strip()
        output_path = self.output_edit.text().strip()
        gene_ids = self._parse_gene_ids()
        promoter_length = self.length_spin.value()

        errors = []
        if not genome_path or not os.path.isfile(genome_path):
            errors.append("Genome FASTA file is missing or does not exist.")
        if not gff_path or not os.path.isfile(gff_path):
            errors.append("GFF3 annotation file is missing or does not exist.")
        if not gene_ids:
            errors.append("Please enter at least one gene ID.")
        if not output_path:
            errors.append("Please choose an output FASTA path.")

        if errors:
            QMessageBox.warning(self, "Cannot start extraction", "\n".join(f"• {e}" for e in errors))
            return

        out_dir = os.path.dirname(os.path.abspath(output_path))
        if out_dir and not os.path.isdir(out_dir):
            try:
                os.makedirs(out_dir, exist_ok=True)
            except Exception as e:
                QMessageBox.warning(self, "Cannot create output folder", str(e))
                return

        self.log_view.clear()
        self.progress_bar.setRange(0, 0)  # busy while loading genome
        self.progress_bar.setValue(0)
        self._set_inputs_enabled(False)
        self.open_folder_btn.setEnabled(False)
        self.status_bar.showMessage("Running extraction...")

        self.worker = ExtractionWorker(
            genome_path, gff_path, gene_ids, output_path, promoter_length,
            strip_version=self.strip_version_check.isChecked()
        )
        self.worker.log.connect(self._append_log)
        self.worker.progress.connect(self._on_progress)
        self.worker.finished_ok.connect(self._on_finished)
        self.worker.failed.connect(self._on_failed)
        self.worker.start()

    def _on_progress(self, current, total):
        if total <= 0:
            self.progress_bar.setRange(0, 0)  # indeterminate
        else:
            self.progress_bar.setRange(0, total)
            self.progress_bar.setValue(current)

    def _cancel_extraction(self):
        if self.worker:
            self.worker.cancel()
        self.cancel_btn.setEnabled(False)
        self.status_bar.showMessage("Cancelling...")

    def _on_finished(self, count, missing, output_path):
        self._set_inputs_enabled(True)
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(100)
        self.open_folder_btn.setEnabled(True)
        self.status_bar.showMessage(f"Done — {count} promoter(s) extracted.")

        self._append_log(f"\nExtracted {count} promoter(s) -> {output_path}")
        if missing:
            self._append_log(f"\n{len(missing)} gene(s) not found:")
            for g in missing:
                self._append_log(f"  {g}")

        summary = f"Extracted {count} promoter sequence(s).\nSaved to:\n{output_path}"
        if missing:
            summary += f"\n\n{len(missing)} gene ID(s) were not found (see log for details)."
        QMessageBox.information(self, "Extraction complete", summary)

    def _on_failed(self, error_text):
        self._set_inputs_enabled(True)
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        self.status_bar.showMessage("Failed.")
        self._append_log("\nERROR:\n" + error_text)
        QMessageBox.critical(self, "Extraction failed", error_text.strip().splitlines()[-1])

    def closeEvent(self, event):
        if self.worker and self.worker.isRunning():
            self.worker.cancel()
            self.worker.wait(2000)
        event.accept()


def main():
    app = QApplication(sys.argv)
    app.setApplicationName(APP_TITLE)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()