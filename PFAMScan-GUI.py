#!/usr/bin/env python3
import subprocess
import sys
import os
import pandas as pd
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QLineEdit, QPushButton, QTextEdit, QGroupBox,
    QFileDialog, QSpinBox, QStatusBar, QTabWidget, QProgressBar,
    QMessageBox, QSplitter, QMenuBar, QCheckBox
)
from PySide6.QtCore import Qt, QThread, Signal, QObject
from PySide6.QtGui import QFont, QTextCursor, QPalette, QColor, QAction, QActionGroup

# ====================== WORKER THREAD ======================
class PfamWorker(QObject):
    log_signal = Signal(str)
    progress_signal = Signal(int)
    finished_signal = Signal()
    error_signal = Signal(str)
    file_output_signal = Signal(str, str)

    def __init__(self, config):
        super().__init__()
        self.config = config

    def run(self):
        try:
            self.log_signal.emit("Starting Pfam Scan...")
            self.log_signal.emit(f"FASTA: {self.config['fasta_file']}")
            self.log_signal.emit(f"Pfam DB: {self.config['pfam_dir']}")
            self.log_signal.emit(f"Output Dir: {self.config['output_dir']}")

            if not os.path.exists(self.config['fasta_file']):
                raise FileNotFoundError(f"FASTA file not found: {self.config['fasta_file']}")
            if not os.path.exists(self.config['pfam_dir']):
                raise FileNotFoundError(f"Pfam directory not found: {self.config['pfam_dir']}")

            os.makedirs(self.config['output_dir'], exist_ok=True)
            self.check_pfam_installed()
            self.progress_signal.emit(20)

            outputs = self.run_pfam_scan()
            self.progress_signal.emit(80)

            parsed_csv = os.path.join(self.config['output_dir'], f"{self.config['output_base']}.parsed.csv")
            self.parse_pfam_output(outputs['stdout'], parsed_csv)
            self.file_output_signal.emit("pfam_output", parsed_csv)
            self.progress_signal.emit(100)

            self.log_signal.emit("\n[✓] Pfam Scan completed successfully!")

        except Exception as e:
            self.error_signal.emit(str(e))
        finally:
            self.finished_signal.emit()

    def check_pfam_installed(self):
        try:
            subprocess.run(["pfam_scan.pl", "-h"], check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            self.log_signal.emit("[✓] pfam_scan.pl detected.")
        except Exception:
            raise RuntimeError("pfam_scan.pl not found. Install PfamScan or add to PATH.")

    def run_pfam_scan(self):
        output_file = os.path.join(self.config['output_dir'], f"{self.config['output_base']}.pfam.txt")

        cmd = [
            "pfam_scan.pl",
            "-fasta", self.config['fasta_file'],
            "-dir", self.config['pfam_dir'],
            "-outfile", output_file,
            "-cpu", str(self.config['cpu']),
            "-e_dom", str(self.config['dom_evalue']),
            "-e_seq", str(self.config['seq_evalue'])
        ]

        # ✅ Add --as if user selected Predict Active Sites
        if self.config.get("predict_as", False):
            cmd.append("-as")

        self.log_signal.emit(f"[i] Running PfamScan:\n{' '.join(cmd)}\n")

        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for line in process.stdout:
            self.log_signal.emit(line.strip())
        process.wait()

        if process.returncode != 0:
            raise RuntimeError("PfamScan failed to execute properly.")

        self.log_signal.emit(f"\n[✓] PfamScan output: {output_file}")
        return {"stdout": output_file}

    def parse_pfam_output(self, output_file, out_csv):
        self.log_signal.emit(f"[i] Parsing PfamScan output: {output_file}")
        cols = [
            "seq_id", "alignment_start", "alignment_end", "envelope_start",
            "envelope_end", "hmm_acc", "hmm_name", "type", "hmm_start",
            "hmm_end", "hmm_length", "bit_score", "E-value",
            "significance", "clan", "predicted_active_site_residues"
        ]
        df = pd.read_csv(output_file, sep=r"\s+", names=cols, comment="#", engine="python")
        df.to_csv(out_csv, index=False)
        self.log_signal.emit(f"[✓] Parsed PfamScan results saved to: {out_csv}")


# ====================== MAIN GUI ======================
class PfamGUI(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("PFAM Scanner")
        self.setGeometry(100, 100, 900, 700)
        self.create_menu_bar()

        central = QWidget()
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)
        tabs = QTabWidget()
        layout.addWidget(tabs)

        # Input Tab
        input_tab = QWidget()
        tabs.addTab(input_tab, "Input")
        input_layout = QVBoxLayout(input_tab)

        group = QGroupBox("Pfam Configuration")
        form_layout = QVBoxLayout(group)
        input_layout.addWidget(group)

        # FASTA
        fasta_layout = QHBoxLayout()
        self.fasta_input = QLineEdit()
        self.fasta_input.setPlaceholderText("/path/to/query.fasta")
        fasta_btn = QPushButton("Browse")
        fasta_btn.clicked.connect(self.browse_fasta)
        fasta_layout.addWidget(QLabel("FASTA File:"))
        fasta_layout.addWidget(self.fasta_input)
        fasta_layout.addWidget(fasta_btn)
        form_layout.addLayout(fasta_layout)

        # Pfam dir
        pfam_layout = QHBoxLayout()
        self.pfam_dir_input = QLineEdit()
        self.pfam_dir_input.setPlaceholderText("/path/to/PfamDB")
        pfam_btn = QPushButton("Browse")
        pfam_btn.clicked.connect(self.browse_pfam_dir)
        pfam_layout.addWidget(QLabel("Pfam DB Directory:"))
        pfam_layout.addWidget(self.pfam_dir_input)
        pfam_layout.addWidget(pfam_btn)
        form_layout.addLayout(pfam_layout)

        # Output dir
        out_layout = QHBoxLayout()
        self.output_input = QLineEdit(os.getcwd())
        out_btn = QPushButton("Browse")
        out_btn.clicked.connect(self.browse_output_dir)
        out_layout.addWidget(QLabel("Output Directory:"))
        out_layout.addWidget(self.output_input)
        out_layout.addWidget(out_btn)
        form_layout.addLayout(out_layout)

        # Output name
        base_layout = QHBoxLayout()
        self.base_input = QLineEdit("pfamscan_results")
        base_layout.addWidget(QLabel("Output Base Name:"))
        base_layout.addWidget(self.base_input)
        form_layout.addLayout(base_layout)

        # Parameters
        param_layout = QHBoxLayout()
        self.seq_eval = QLineEdit("1e-5")
        self.dom_eval = QLineEdit("1e-5")
        self.cpu_spin = QSpinBox()
        self.cpu_spin.setRange(1, os.cpu_count())
        self.cpu_spin.setValue(4)
        param_layout.addWidget(QLabel("Seq E-value:"))
        param_layout.addWidget(self.seq_eval)
        param_layout.addWidget(QLabel("Dom E-value:"))
        param_layout.addWidget(self.dom_eval)
        param_layout.addWidget(QLabel("CPU:"))
        param_layout.addWidget(self.cpu_spin)
        form_layout.addLayout(param_layout)

        # ✅ Active Site Checkbox
        self.as_checkbox = QCheckBox("Predict Active Sites (-as)")
        form_layout.addWidget(self.as_checkbox)

        # Run button
        run_btn = QPushButton("Run Pfam Scan")
        run_btn.clicked.connect(self.run_scan)
        form_layout.addWidget(run_btn)

        # Output Tab
        output_tab = QWidget()
        tabs.addTab(output_tab, "Output")
        output_layout = QVBoxLayout(output_tab)

        self.console = QTextEdit()
        self.console.setReadOnly(True)
        self.console.setFont(QFont("Courier", 10))

        self.progress = QProgressBar()
        output_layout.addWidget(self.progress)
        output_layout.addWidget(self.console)

        # Status bar
        self.status = QStatusBar()
        self.setStatusBar(self.status)
        self.worker = None
        self.thread = None
        self.set_dark_theme()

    def create_menu_bar(self):
        menubar = QMenuBar()
        self.setMenuBar(menubar)
        help_menu = menubar.addMenu("Help")
        about = QAction("About", self)
        about.triggered.connect(lambda: QMessageBox.about(
            self, "About", "PfamScan GUI v1.1\nSupports Active Site Prediction (-as option)."))
        help_menu.addAction(about)

    def set_dark_theme(self):
        palette = QPalette()
        palette.setColor(QPalette.Window, QColor(45, 45, 45))
        palette.setColor(QPalette.WindowText, Qt.white)
        palette.setColor(QPalette.Base, QColor(25, 25, 25))
        palette.setColor(QPalette.Text, Qt.white)
        QApplication.instance().setPalette(palette)

    def browse_fasta(self):
        f, _ = QFileDialog.getOpenFileName(self, "Select FASTA", "", "FASTA (*.fasta *.fa *.faa)")
        if f:
            self.fasta_input.setText(f)

    def browse_pfam_dir(self):
        d = QFileDialog.getExistingDirectory(self, "Select Pfam Directory")
        if d:
            self.pfam_dir_input.setText(d)

    def browse_output_dir(self):
        d = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if d:
            self.output_input.setText(d)

    def run_scan(self):
        config = {
            'fasta_file': self.fasta_input.text(),
            'pfam_dir': self.pfam_dir_input.text(),
            'output_dir': self.output_input.text(),
            'output_base': self.base_input.text(),
            'seq_evalue': self.seq_eval.text(),
            'dom_evalue': self.dom_eval.text(),
            'cpu': self.cpu_spin.value(),
            'predict_as': self.as_checkbox.isChecked()   # ✅ include AS flag
        }

        self.worker = PfamWorker(config)
        self.thread = QThread()
        self.worker.moveToThread(self.thread)

        self.worker.log_signal.connect(self.log)
        self.worker.progress_signal.connect(self.progress.setValue)
        self.worker.finished_signal.connect(self.finish)
        self.worker.error_signal.connect(self.error)

        self.thread.started.connect(self.worker.run)
        self.thread.start()
        self.status.showMessage("Running Pfam Scan...")

    def log(self, msg):
        self.console.append(msg)
        self.console.moveCursor(QTextCursor.End)

    def finish(self):
        self.status.showMessage("Pfam Scan Completed.")
        QMessageBox.information(self, "Done", "Pfam Scan completed successfully!")

    def error(self, msg):
        self.console.append(f"\n[ERROR] {msg}")
        QMessageBox.critical(self, "Error", msg)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = PfamGUI()
    window.show()
    sys.exit(app.exec())
