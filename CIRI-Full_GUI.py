#!/usr/bin/env python3

import os
import re
import sys
import subprocess
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, 
    QGridLayout, QLabel, QLineEdit, QPushButton, QFileDialog, 
    QSpinBox, QGroupBox, QTextEdit, QMessageBox, QProgressBar
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal
from PyQt6.QtGui import QFont, QAction

# ---------------- WORKER THREAD ---------------- #
class PipelineWorker(QThread):
    log_signal = pyqtSignal(str)
    progress_signal = pyqtSignal(int)
    finished_signal = pyqtSignal(bool, str)

    def __init__(self, fastq_dir, ref, gtf, outdir, threads):
        super().__init__()
        self.fastq_dir = fastq_dir
        self.ref = ref
        self.gtf = gtf
        self.outdir = outdir
        self.threads = threads
        self._is_running = True
        self.current_process = None

    def run(self):
        try:
            os.makedirs(self.outdir, exist_ok=True)

            fastq_files = [f for f in os.listdir(self.fastq_dir) if f.endswith(".fastq.gz")]
            if not fastq_files:
                self.finished_signal.emit(False, "No .fastq.gz files found in the selected directory.")
                return

            pairs = {}
            for f in fastq_files:
                m = re.match(r"(.+?)_(1|2)\.fastq\.gz$", f)
                if not m:
                    self.log_signal.emit(f"[WARNING] Skipping unrecognized file: {f}")
                    continue
                sample, read = m.groups()
                pairs.setdefault(sample, {})[read] = f

            valid_samples = {s: r for s, r in pairs.items() if "1" in r and "2" in r}
            total_samples = len(valid_samples)

            if total_samples == 0:
                self.finished_signal.emit(False, "No valid paired fastq files found.")
                return

            self.log_signal.emit(f"[INFO] Found {total_samples} valid paired samples. Starting pipeline...\n")

            for index, (sample, reads) in enumerate(sorted(valid_samples.items())):
                if not self._is_running:
                    self.log_signal.emit("\n[INFO] Pipeline stopped by user.")
                    break

                r1 = reads["1"]
                r2 = reads["2"]
                sample_outdir = os.path.join(self.outdir, sample)

                if os.path.exists(sample_outdir):
                    self.log_signal.emit(f"[INFO] Skipping {sample}: output directory already exists.")
                    self.progress_signal.emit(int(((index + 1) / total_samples) * 100))
                    continue

                cmd = [
                    "CIRI-full", "Pipeline",
                    "-1", os.path.join(self.fastq_dir, r1),
                    "-2", os.path.join(self.fastq_dir, r2),
                    "-r", self.ref,
                    "-d", sample_outdir,
                    "-o", sample,
                    "-t", str(self.threads),
                    "-0"
                ]

                if self.gtf:
                    cmd.extend(["-a", self.gtf])

                self.log_signal.emit(f"▶ Running CIRI-full for sample: {sample}")
                self.log_signal.emit(f"  Command: {' '.join(cmd)}\n")

                self.current_process = subprocess.Popen(
                    cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1
                )

                for line in self.current_process.stdout:
                    self.log_signal.emit(line.strip())
                
                self.current_process.wait()

                if self.current_process.returncode != 0 and self._is_running:
                    self.log_signal.emit(f"\n[ERROR] CIRI-full failed for sample {sample} with exit code {self.current_process.returncode}")
                
                self.progress_signal.emit(int(((index + 1) / total_samples) * 100))

            if self._is_running:
                self.finished_signal.emit(True, "All samples processed successfully.")

        except Exception as e:
            self.finished_signal.emit(False, f"An unexpected error occurred: {str(e)}")

    def stop(self):
        self._is_running = False
        if self.current_process:
            self.current_process.terminate()


# ---------------- MAIN GUI WINDOW ---------------- #
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("CIRI-full Batch Pipeline")
        self.resize(800, 650)
        self.worker = None

        self._create_menu_bar()

        # Main Widget and Layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QVBoxLayout(central_widget)

        # 1. Inputs Group
        input_group = QGroupBox("Configuration")
        input_layout = QGridLayout()
        input_group.setLayout(input_layout)

        # FASTQ Dir
        self.fastq_input = QLineEdit()
        self.fastq_btn = QPushButton("Browse")
        self.fastq_btn.clicked.connect(self.browse_fastq_dir)
        input_layout.addWidget(QLabel("Cleaned FASTQ Directory:"), 0, 0)
        input_layout.addWidget(self.fastq_input, 0, 1)
        input_layout.addWidget(self.fastq_btn, 0, 2)

        # Reference FASTA
        self.ref_input = QLineEdit()
        self.ref_btn = QPushButton("Browse")
        self.ref_btn.clicked.connect(self.browse_ref_file)
        input_layout.addWidget(QLabel("Reference FASTA:"), 1, 0)
        input_layout.addWidget(self.ref_input, 1, 1)
        input_layout.addWidget(self.ref_btn, 1, 2)

        # Annotation GTF
        self.gtf_input = QLineEdit()
        self.gtf_input.setPlaceholderText("Optional but recommended")
        self.gtf_btn = QPushButton("Browse")
        self.gtf_btn.clicked.connect(self.browse_gtf_file)
        input_layout.addWidget(QLabel("Annotation GTF:"), 2, 0)
        input_layout.addWidget(self.gtf_input, 2, 1)
        input_layout.addWidget(self.gtf_btn, 2, 2)

        # Output Dir
        self.out_input = QLineEdit()
        self.out_btn = QPushButton("Browse")
        self.out_btn.clicked.connect(self.browse_out_dir)
        input_layout.addWidget(QLabel("Output Directory:"), 3, 0)
        input_layout.addWidget(self.out_input, 3, 1)
        input_layout.addWidget(self.out_btn, 3, 2)

        # Threads
        self.thread_spin = QSpinBox()
        self.thread_spin.setRange(1, 128)
        self.thread_spin.setValue(10)
        input_layout.addWidget(QLabel("Threads:"), 4, 0)
        input_layout.addWidget(self.thread_spin, 4, 1)

        main_layout.addWidget(input_group)

        # 2. Controls Group
        control_layout = QHBoxLayout()
        self.run_btn = QPushButton("Run Pipeline")
        self.run_btn.setStyleSheet("background-color: #2E8B57; color: white; font-weight: bold; padding: 8px;")
        self.run_btn.clicked.connect(self.start_pipeline)
        
        self.stop_btn = QPushButton("Stop")
        self.stop_btn.setStyleSheet("background-color: #B22222; color: white; font-weight: bold; padding: 8px;")
        self.stop_btn.clicked.connect(self.stop_pipeline)
        self.stop_btn.setEnabled(False)

        control_layout.addWidget(self.run_btn)
        control_layout.addWidget(self.stop_btn)
        main_layout.addLayout(control_layout)

        # 3. Progress Bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setValue(0)
        main_layout.addWidget(self.progress_bar)

        # 4. Console/Log Output
        log_group = QGroupBox("Pipeline Log")
        log_layout = QVBoxLayout()
        log_group.setLayout(log_layout)
        self.log_console = QTextEdit()
        self.log_console.setReadOnly(True)
        self.log_console.setFont(QFont("Courier", 10))
        self.log_console.setStyleSheet("background-color: #1E1E1E; color: #DCDCAA;")
        log_layout.addWidget(self.log_console)
        
        main_layout.addWidget(log_group)

    # --- Menu Bar Setup ---
    def _create_menu_bar(self):
        menu_bar = self.menuBar()

        # File Menu
        file_menu = menu_bar.addMenu("File")
        
        clear_log_action = QAction("Clear Log", self)
        clear_log_action.setShortcut("Ctrl+L")
        clear_log_action.triggered.connect(self.clear_log)
        file_menu.addAction(clear_log_action)
        
        file_menu.addSeparator()
        
        exit_action = QAction("Exit", self)
        exit_action.setShortcut("Ctrl+Q")
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)

        # Help Menu
        help_menu = menu_bar.addMenu("Help")
        
        about_action = QAction("About", self)
        about_action.triggered.connect(self.show_about)
        help_menu.addAction(about_action)

    # --- Application Logic ---
    def clear_log(self):
        self.log_console.clear()

    def show_about(self):
        QMessageBox.about(
            self, 
            "About CIRI-full Pipeline", 
            "<h3>CIRI-full GUI Wrapper</h3>"
            "<p>A graphical interface for running the CIRI-full circRNA pipeline in batch mode.</p>"
            "<p>Features include background multithreading, real-time logging, and graceful process termination.</p>"
        )

    # --- Directory / File Browsers ---
    def browse_fastq_dir(self):
        d = QFileDialog.getExistingDirectory(self, "Select FASTQ Directory")
        if d: self.fastq_input.setText(d)

    def browse_ref_file(self):
        f, _ = QFileDialog.getOpenFileName(self, "Select Reference FASTA", "", "FASTA Files (*.fa *.fasta);;All Files (*)")
        if f: self.ref_input.setText(f)

    def browse_gtf_file(self):
        f, _ = QFileDialog.getOpenFileName(self, "Select Annotation GTF", "", "GTF Files (*.gtf);;All Files (*)")
        if f: self.gtf_input.setText(f)

    def browse_out_dir(self):
        d = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if d: self.out_input.setText(d)

    # --- Pipeline Logic ---
    def append_log(self, text):
        self.log_console.append(text)
        scrollbar = self.log_console.verticalScrollBar()
        scrollbar.setValue(scrollbar.maximum())

    def start_pipeline(self):
        fastq = self.fastq_input.text().strip()
        ref = self.ref_input.text().strip()
        gtf = self.gtf_input.text().strip()
        out = self.out_input.text().strip()
        threads = self.thread_spin.value()

        if not fastq or not ref or not out:
            QMessageBox.warning(self, "Missing Inputs", "Please provide FASTQ Directory, Reference FASTA, and Output Directory.")
            return

        self.run_btn.setEnabled(False)
        self.stop_btn.setEnabled(True)
        self.progress_bar.setValue(0)
        self.log_console.clear()
        self.append_log("[SYSTEM] Starting pipeline...\n")

        self.worker = PipelineWorker(fastq, ref, gtf, out, threads)
        self.worker.log_signal.connect(self.append_log)
        self.worker.progress_signal.connect(self.progress_bar.setValue)
        self.worker.finished_signal.connect(self.pipeline_finished)
        self.worker.start()

    def stop_pipeline(self):
        if self.worker and self.worker.isRunning():
            self.append_log("\n[SYSTEM] Sending stop signal... Waiting for current process to terminate.")
            self.worker.stop()
            self.stop_btn.setEnabled(False)

    def pipeline_finished(self, success, message):
        self.run_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)
        
        if success:
            self.append_log(f"\n[SUCCESS] {message}")
            QMessageBox.information(self, "Complete", message)
        else:
            self.append_log(f"\n[FINISHED/FAILED] {message}")
            QMessageBox.warning(self, "Warning", message)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    window = MainWindow()
    window.show()
    sys.exit(app.exec())