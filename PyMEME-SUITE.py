import sys
import os
import subprocess
from PyQt6.QtWidgets import (
    QApplication, QWidget, QTabWidget, QVBoxLayout, QHBoxLayout, QLabel,
    QPushButton, QFileDialog, QLineEdit, QTextEdit, QSpinBox, QComboBox, QMessageBox
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal


class MemeRunner(QThread):
    """Threaded process runner for MEME/TOMTOM commands"""
    output_signal = pyqtSignal(str)
    finished_signal = pyqtSignal(bool)

    def __init__(self, cmd):
        super().__init__()
        self.cmd = cmd

    def run(self):
        process = subprocess.Popen(
            self.cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True
        )
        for line in iter(process.stdout.readline, ""):
            if not line:
                break
            self.output_signal.emit(line)
        process.wait()
        self.finished_signal.emit(process.returncode == 0)


class MemeSuiteApp(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("MEME Suite GUI - by Joydeep Shaw 🧬")
        self.resize(950, 600)
        layout = QVBoxLayout()
        self.tabs = QTabWidget()

        self.meme_tab = self.build_meme_tab()
        self.tomtom_tab = self.build_tomtom_tab()

        self.tabs.addTab(self.meme_tab, "Run MEME")
        self.tabs.addTab(self.tomtom_tab, "Run TOMTOM")
        layout.addWidget(self.tabs)
        self.setLayout(layout)

    # ------------------------ MEME TAB ------------------------
    def build_meme_tab(self):
        widget = QWidget()
        layout = QVBoxLayout()

        # FASTA input
        fasta_layout = QHBoxLayout()
        self.fasta_path = QLineEdit()
        fasta_btn = QPushButton("Select FASTA File")
        fasta_btn.clicked.connect(self.select_fasta)
        fasta_layout.addWidget(QLabel("Input FASTA:"))
        fasta_layout.addWidget(self.fasta_path)
        fasta_layout.addWidget(fasta_btn)
        layout.addLayout(fasta_layout)

        # Output directory
        out_layout = QHBoxLayout()
        self.out_path = QLineEdit("meme_out")
        out_btn = QPushButton("Select Output Dir")
        out_btn.clicked.connect(self.select_output)
        out_layout.addWidget(QLabel("Output Dir:"))
        out_layout.addWidget(self.out_path)
        out_layout.addWidget(out_btn)
        layout.addLayout(out_layout)

        # Parameters section
        param_layout = QHBoxLayout()

        # Sequence type dropdown
        self.seq_type_box = QComboBox()
        self.seq_type_box.addItems(["DNA", "RNA", "Protein"])

        # Model dropdown
        self.model_box = QComboBox()
        self.model_box.addItems(["zoops", "oops", "anr"])

        # Motif numbers
        self.nmotifs_spin = QSpinBox()
        self.nmotifs_spin.setRange(1, 100)
        self.nmotifs_spin.setValue(5)

        # Widths
        self.minw_spin = QSpinBox()
        self.minw_spin.setRange(4, 50)
        self.minw_spin.setValue(6)

        self.maxw_spin = QSpinBox()
        self.maxw_spin.setRange(6, 100)
        self.maxw_spin.setValue(15)

        # Threads
        self.threads_spin = QSpinBox()
        self.threads_spin.setRange(1, 32)
        self.threads_spin.setValue(4)

        # Add widgets
        param_layout.addWidget(QLabel("Seq Type:"))
        param_layout.addWidget(self.seq_type_box)
        param_layout.addWidget(QLabel("Model:"))
        param_layout.addWidget(self.model_box)
        param_layout.addWidget(QLabel("Motifs:"))
        param_layout.addWidget(self.nmotifs_spin)
        param_layout.addWidget(QLabel("Min W:"))
        param_layout.addWidget(self.minw_spin)
        param_layout.addWidget(QLabel("Max W:"))
        param_layout.addWidget(self.maxw_spin)
        param_layout.addWidget(QLabel("Threads:"))
        param_layout.addWidget(self.threads_spin)

        layout.addLayout(param_layout)

        # Run button
        run_btn = QPushButton("▶ Run MEME")
        run_btn.clicked.connect(self.run_meme)
        layout.addWidget(run_btn)

        # Output console
        self.console = QTextEdit()
        self.console.setReadOnly(True)
        layout.addWidget(QLabel("Execution Log:"))
        layout.addWidget(self.console)

        widget.setLayout(layout)
        return widget

    # ------------------------ TOMTOM TAB ------------------------
    def build_tomtom_tab(self):
        widget = QWidget()
        layout = QVBoxLayout()

        meme_layout = QHBoxLayout()
        self.meme_file = QLineEdit()
        meme_btn = QPushButton("Select MEME File")
        meme_btn.clicked.connect(self.select_meme)
        meme_layout.addWidget(QLabel("MEME Motif File:"))
        meme_layout.addWidget(self.meme_file)
        meme_layout.addWidget(meme_btn)
        layout.addLayout(meme_layout)

        db_layout = QHBoxLayout()
        self.db_file = QLineEdit()
        db_btn = QPushButton("Select Motif Database")
        db_btn.clicked.connect(self.select_db)
        db_layout.addWidget(QLabel("Motif Database:"))
        db_layout.addWidget(self.db_file)
        db_layout.addWidget(db_btn)
        layout.addLayout(db_layout)

        self.tomtom_out = QLineEdit("tomtom_out")
        layout.addWidget(QLabel("Output Directory:"))
        layout.addWidget(self.tomtom_out)

        run_tomtom_btn = QPushButton("▶ Run TOMTOM")
        run_tomtom_btn.clicked.connect(self.run_tomtom)
        layout.addWidget(run_tomtom_btn)

        self.tomtom_console = QTextEdit()
        self.tomtom_console.setReadOnly(True)
        layout.addWidget(QLabel("Execution Log:"))
        layout.addWidget(self.tomtom_console)

        widget.setLayout(layout)
        return widget

    # ------------------------ File Selectors ------------------------
    def select_fasta(self):
        file, _ = QFileDialog.getOpenFileName(self, "Select FASTA File", "", "FASTA files (*.fa *.fasta)")
        if file:
            self.fasta_path.setText(file)

    def select_output(self):
        dir_ = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if dir_:
            self.out_path.setText(dir_)

    def select_meme(self):
        file, _ = QFileDialog.getOpenFileName(self, "Select MEME Motif File", "", "MEME files (*.txt *.meme)")
        if file:
            self.meme_file.setText(file)

    def select_db(self):
        file, _ = QFileDialog.getOpenFileName(self, "Select Motif Database", "", "MEME motif db (*.meme)")
        if file:
            self.db_file.setText(file)

    # ------------------------ Command Execution ------------------------
    def run_meme(self):
        fasta = self.fasta_path.text().strip()
        outdir = self.out_path.text().strip()
        if not os.path.exists(fasta):
            QMessageBox.critical(self, "Error", "Please select a valid FASTA file.")
            return

        seq_type = self.seq_type_box.currentText().lower()
        cmd = (
            f"meme {fasta} -oc {outdir} "
            f"-{seq_type} "
            f"-mod {self.model_box.currentText()} "
            f"-nmotifs {self.nmotifs_spin.value()} "
            f"-minw {self.minw_spin.value()} "
            f"-maxw {self.maxw_spin.value()} "
            f"-p {self.threads_spin.value()}"
        )

        self.console.append(f"\n🔹 Running command:\n{cmd}\n")
        self.thread = MemeRunner(cmd)
        self.thread.output_signal.connect(self.console.append)
        self.thread.finished_signal.connect(lambda ok: self.on_finished(ok, outdir))
        self.thread.start()

    def run_tomtom(self):
        meme_file = self.meme_file.text().strip()
        db_file = self.db_file.text().strip()
        outdir = self.tomtom_out.text().strip()

        if not (os.path.exists(meme_file) and os.path.exists(db_file)):
            QMessageBox.critical(self, "Error", "Please select valid MEME and database files.")
            return

        cmd = f"tomtom -oc {outdir} {meme_file} {db_file}"
        self.tomtom_console.append(f"\n🔹 Running command:\n{cmd}\n")
        self.thread = MemeRunner(cmd)
        self.thread.output_signal.connect(self.tomtom_console.append)
        self.thread.finished_signal.connect(lambda ok: self.on_finished(ok, outdir))
        self.thread.start()

    def on_finished(self, ok, outdir):
        if ok:
            QMessageBox.information(self, "✅ Success", f"Process completed successfully!\nResults saved in:\n{outdir}")
        else:
            QMessageBox.warning(self, "⚠️ Failed", "Process encountered an error.")


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MemeSuiteApp()
    window.show()
    sys.exit(app.exec())
