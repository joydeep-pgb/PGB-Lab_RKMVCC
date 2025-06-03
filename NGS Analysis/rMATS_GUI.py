import os
import sys
import subprocess
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QLineEdit, QPushButton, QGroupBox, QSpinBox, QCheckBox,
    QFileDialog, QTextEdit, QComboBox, QProgressBar
)
from PySide6.QtCore import QProcess, Qt
from PySide6.QtGui import QColor, QPalette

class rMATSGUI(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("rMATS Alternative Splicing Analyzer")
        self.setGeometry(100, 100, 800, 600)

        # Central Widget and Layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QVBoxLayout(central_widget)

        # Input Parameters Group
        input_group = QGroupBox("Input Parameters")
        input_layout = QVBoxLayout(input_group)

        # Batch 1 (always required)
        b1_layout = QHBoxLayout()
        self.b1_label = QLabel("Batch 1 File*:")
        b1_layout.addWidget(self.b1_label)
        self.b1_edit = QLineEdit()
        b1_layout.addWidget(self.b1_edit)
        self.b1_button = QPushButton("Browse")
        self.b1_button.clicked.connect(lambda: self.browse_file(self.b1_edit))
        b1_layout.addWidget(self.b1_button)
        input_layout.addLayout(b1_layout)

        # Batch 2 (optional when statoff enabled)
        b2_layout = QHBoxLayout()
        self.b2_label = QLabel("Batch 2 File*:")
        b2_layout.addWidget(self.b2_label)
        self.b2_edit = QLineEdit()
        b2_layout.addWidget(self.b2_edit)
        self.b2_button = QPushButton("Browse")
        self.b2_button.clicked.connect(lambda: self.browse_file(self.b2_edit))
        b2_layout.addWidget(self.b2_button)
        input_layout.addLayout(b2_layout)

        # GTF File (always required)
        gtf_layout = QHBoxLayout()
        self.gtf_label = QLabel("GTF File*:")
        gtf_layout.addWidget(self.gtf_label)
        self.gtf_edit = QLineEdit()
        gtf_layout.addWidget(self.gtf_edit)
        self.gtf_button = QPushButton("Browse")
        self.gtf_button.clicked.connect(lambda: self.browse_file(self.gtf_edit))
        gtf_layout.addWidget(self.gtf_button)
        input_layout.addLayout(gtf_layout)

        # Output Directory (always required)
        out_layout = QHBoxLayout()
        self.out_label = QLabel("Output Directory*:")
        out_layout.addWidget(self.out_label)
        self.out_edit = QLineEdit()
        out_layout.addWidget(self.out_edit)
        self.out_button = QPushButton("Browse")
        self.out_button.clicked.connect(self.browse_directory)
        out_layout.addWidget(self.out_button)
        input_layout.addLayout(out_layout)

        # Temp Directory (always required)
        temp_layout = QHBoxLayout()
        self.temp_label = QLabel("Temp Directory*:")
        temp_layout.addWidget(self.temp_label)
        self.temp_edit = QLineEdit()
        temp_layout.addWidget(self.temp_edit)
        self.temp_button = QPushButton("Browse")
        self.temp_button.clicked.connect(self.browse_temp_directory)
        temp_layout.addWidget(self.temp_button)
        input_layout.addLayout(temp_layout)

        # Options Group
        options_group = QGroupBox("Analysis Options")
        options_layout = QVBoxLayout(options_group)

        # Read Type
        type_layout = QHBoxLayout()
        type_layout.addWidget(QLabel("Read Type*:"))
        self.type_combo = QComboBox()
        self.type_combo.addItems(["paired", "single"])
        self.type_combo.setCurrentText("paired")
        type_layout.addWidget(self.type_combo)
        options_layout.addLayout(type_layout)

        # Read Length
        read_len_layout = QHBoxLayout()
        read_len_layout.addWidget(QLabel("Read Length*:"))
        self.read_len_spin = QSpinBox()
        self.read_len_spin.setRange(1, 500)
        self.read_len_spin.setValue(295)
        read_len_layout.addWidget(self.read_len_spin)
        options_layout.addLayout(read_len_layout)

        # Threads
        threads_layout = QHBoxLayout()
        threads_layout.addWidget(QLabel("Threads*:"))
        self.threads_spin = QSpinBox()
        self.threads_spin.setRange(1, 64)
        self.threads_spin.setValue(8)
        threads_layout.addWidget(self.threads_spin)
        options_layout.addLayout(threads_layout)

        # Checkboxes
        self.var_read_check = QCheckBox("Variable Read Length")
        self.var_read_check.setChecked(True)
        options_layout.addWidget(self.var_read_check)

        self.statoff_check = QCheckBox("Disable Statistical Analysis (--statoff)")
        self.statoff_check.setChecked(False)  # Disabled by default
        self.statoff_check.stateChanged.connect(self.toggle_b2_requirement)
        options_layout.addWidget(self.statoff_check)

        # Progress and Output
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        main_layout.addWidget(self.progress_bar)

        self.output_text = QTextEdit()
        self.output_text.setReadOnly(True)
        main_layout.addWidget(self.output_text)

        # Run Button
        self.run_button = QPushButton("Run rMATS Analysis")
        self.run_button.clicked.connect(self.run_analysis)
        main_layout.addWidget(self.run_button)

        # Add groups to main layout
        main_layout.addWidget(input_group)
        main_layout.addWidget(options_group)

        # QProcess for running commands
        self.process = QProcess()
        self.process.readyReadStandardOutput.connect(self.handle_stdout)
        self.process.readyReadStandardError.connect(self.handle_stderr)
        self.process.finished.connect(self.analysis_finished)

        # Initial UI update
        self.toggle_b2_requirement()

    def toggle_b2_requirement(self):
        """Update UI based on statoff selection"""
        if self.statoff_check.isChecked():
            self.b2_label.setText("Batch 2 File (optional):")
            # Make background lighter to indicate optional
            pal = self.b2_edit.palette()
            pal.setColor(QPalette.Base, QColor(240, 240, 240))
            self.b2_edit.setPalette(pal)
        else:
            self.b2_label.setText("Batch 2 File*:")
            # Reset to default background
            self.b2_edit.setPalette(QApplication.palette())

    def browse_file(self, line_edit):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select File")
        if file_path:
            line_edit.setText(file_path)

    def browse_directory(self):
        dir_path = QFileDialog.getExistingDirectory(self, "Select Directory")
        if dir_path:
            self.out_edit.setText(dir_path)

    def browse_temp_directory(self):
        dir_path = QFileDialog.getExistingDirectory(self, "Select Temp Directory")
        if dir_path:
            self.temp_edit.setText(dir_path)

    def run_analysis(self):
        # Validate required inputs
        required_fields = [
            ("Batch 1", self.b1_edit.text()),
            ("GTF File", self.gtf_edit.text()),
            ("Output Directory", self.out_edit.text()),
            ("Temp Directory", self.temp_edit.text()),
        ]
        
        # Batch 2 is required only when statoff is disabled
        if not self.statoff_check.isChecked():
            required_fields.append(("Batch 2", self.b2_edit.text()))
        
        # Check for empty required fields
        errors = [name for name, value in required_fields if not value.strip()]
        if errors:
            self.output_text.append(f"ERROR: Missing required fields: {', '.join(errors)}")
            return

        # Build command
        cmd = [
            "rmats.py",
            "--b1", self.b1_edit.text().strip(),
            "--gtf", self.gtf_edit.text().strip(),
            "-t", self.type_combo.currentText(),
            "--readLength", str(self.read_len_spin.value()),
            "--nthread", str(self.threads_spin.value()),
            "--od", self.out_edit.text().strip(),
            "--tmp", self.temp_edit.text().strip()
        ]

        # Add Batch 2 only if provided or required
        if self.b2_edit.text().strip() or not self.statoff_check.isChecked():
            cmd.extend(["--b2", self.b2_edit.text().strip()])

        if self.var_read_check.isChecked():
            cmd.append("--variable-read-length")
        
        if self.statoff_check.isChecked():
            cmd.append("--statoff")

        # Display command
        self.output_text.append("Running command:\n" + " ".join(cmd))
        self.progress_bar.setVisible(True)
        self.run_button.setEnabled(False)

        # Start process
        self.process.start(" ".join(cmd))

    def handle_stdout(self):
        data = self.process.readAllStandardOutput()
        if data:
            self.output_text.append(data.data().decode().strip())

    def handle_stderr(self):
        data = self.process.readAllStandardError()
        if data:
            self.output_text.append(f"<font color='red'>{data.data().decode().strip()}</font>")

    def analysis_finished(self):
        self.progress_bar.setVisible(False)
        self.run_button.setEnabled(True)
        self.output_text.append("\nAnalysis completed!")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = rMATSGUI()
    window.show()
    sys.exit(app.exec())