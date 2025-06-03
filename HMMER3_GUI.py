#!/usr/bin/env python3
import subprocess
import sys
import os
import pandas as pd
import re
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QLineEdit, QPushButton, QTextEdit, QGroupBox,
    QFileDialog, QDoubleSpinBox, QSpinBox, QCheckBox, QStatusBar,
    QTabWidget, QProgressBar, QMessageBox, QSplitter, QFrame,
    QMenuBar
)
from PySide6.QtCore import Qt, QThread, Signal, QObject, QRegularExpression
from PySide6.QtGui import QFont, QTextCursor, QPalette, QColor, QRegularExpressionValidator, QIcon, QKeySequence, QAction, QActionGroup

# ====================== WORKER THREAD ======================
class HmmerWorker(QObject):
    log_signal = Signal(str)
    progress_signal = Signal(int)
    finished_signal = Signal()
    error_signal = Signal(str)
    file_output_signal = Signal(str, str)  # Signal for file paths (file_type, path)

    def __init__(self, config):
        super().__init__()
        self.config = config
        self.is_running = True

    def run(self):
        try:
            self.log_signal.emit("Starting HMMER scan...")
            self.log_signal.emit(f"Profile: {self.config['hmm_profile']}")
            self.log_signal.emit(f"FASTA: {self.config['fasta_file']}")
            self.log_signal.emit(f"Output Directory: {self.config['output_dir']}")
            self.log_signal.emit(f"E-value: {self.config['evalue']}")
            self.log_signal.emit(f"Domain E-value: {self.config['dom_evalue']}")
            
            # Check if files exist
            if not os.path.exists(self.config['hmm_profile']):
                raise FileNotFoundError(f"HMM profile not found: {self.config['hmm_profile']}")
            if not os.path.exists(self.config['fasta_file']):
                raise FileNotFoundError(f"FASTA file not found: {self.config['fasta_file']}")
            
            # Create output directory if it doesn't exist
            os.makedirs(self.config['output_dir'], exist_ok=True)
            
            # Check if HMMER is installed
            self.check_hmmer_installed()
            self.progress_signal.emit(10)
            
            # Press HMM profile
            self.hmmpress_profile()
            self.progress_signal.emit(20)
            
            # Run hmmscan
            output_files = self.run_hmmscan()
            self.progress_signal.emit(70)
            
            # Parse outputs
            tblout_csv = os.path.join(self.config['output_dir'], f"{self.config['output_base']}.tblout.parsed.csv")
            self.parse_tblout(output_files['tblout'], tblout_csv)
            self.file_output_signal.emit("tblout", tblout_csv)
            self.progress_signal.emit(80)
            
            domtblout_csv = os.path.join(self.config['output_dir'], f"{self.config['output_base']}.domtblout.parsed.csv")
            self.parse_domtblout(output_files['domtblout'], domtblout_csv)
            self.file_output_signal.emit("domtblout", domtblout_csv)
            self.progress_signal.emit(90)
            
            self.log_signal.emit("\n[✓] Process completed successfully!")
            self.progress_signal.emit(100)
            
        except Exception as e:
            self.error_signal.emit(str(e))
        finally:
            self.finished_signal.emit()

    def check_hmmer_installed(self):
        try:
            subprocess.run(["hmmscan", "-h"], check=True, 
                           stdout=subprocess.DEVNULL, 
                           stderr=subprocess.DEVNULL)
            self.log_signal.emit("[✓] HMMER is installed")
        except (subprocess.CalledProcessError, FileNotFoundError):
            raise RuntimeError("HMMER not found. Install with conda: `conda install -c bioconda hmmer`")

    def hmmpress_profile(self):
        hmm_file = self.config['hmm_profile']
        required_exts = ['.h3m', '.h3i', '.h3f', '.h3p']
        if not all(os.path.exists(hmm_file + ext) for ext in required_exts):
            self.log_signal.emit(f"[i] Pressing HMM profile: {hmm_file}")
            process = subprocess.Popen(["hmmpress", hmm_file], 
                                      stdout=subprocess.PIPE, 
                                      stderr=subprocess.STDOUT,
                                      text=True)
            
            for line in process.stdout:
                self.log_signal.emit(line.strip())
                
            process.wait()
            if process.returncode != 0:
                raise RuntimeError("hmmpress failed to execute")
            self.log_signal.emit("[✓] HMM profile pressed successfully")
        else:
            self.log_signal.emit("[✓] HMM profile already pressed")

    def run_hmmscan(self):
        outputs = {
            'domtblout': os.path.join(self.config['output_dir'], f"{self.config['output_base']}.domtblout"),
            'tblout': os.path.join(self.config['output_dir'], f"{self.config['output_base']}.tblout"),
            'stdout': os.path.join(self.config['output_dir'], f"{self.config['output_base']}.full.txt")
        }

        cmd = [
            "hmmscan",
            "--cpu", str(self.config['cpu']),
            "-E", str(self.config['evalue']),
            "--domE", str(self.config['dom_evalue']),
            "--tblout", outputs['tblout'],
            "--domtblout", outputs['domtblout']
        ] + self.config['extra_options'] + [self.config['hmm_profile'], self.config['fasta_file']]

        self.log_signal.emit(f"[i] Running hmmscan:\n{' '.join(cmd)}\n")
        
        with open(outputs['stdout'], 'w') as outfile:
            process = subprocess.Popen(cmd, 
                                      stdout=subprocess.PIPE, 
                                      stderr=subprocess.STDOUT,
                                      text=True,
                                      bufsize=1)
            
            for line in process.stdout:
                outfile.write(line)
                self.log_signal.emit(line.strip())
                
            process.wait()
            if process.returncode != 0:
                raise RuntimeError("hmmscan failed to execute")
        
        self.log_signal.emit("\n[✓] Scan complete. Output files:")
        for key, path in outputs.items():
            self.log_signal.emit(f"   - {key}: {path}")
            self.file_output_signal.emit(key, path)
        
        return outputs

    def parse_tblout(self, tbl_file, out_csv):
        self.log_signal.emit(f"[i] Parsing tblout: {tbl_file}")
        records = []
        with open(tbl_file) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split()
                if len(fields) >= 18:
                    records.append({
                        "target_name": fields[0],
                        "target_accession": fields[1],
                        "query_name": fields[2],
                        "query_accession": fields[3],
                        "E-value": float(fields[4]),
                        "score": float(fields[5]),
                        "bias": float(fields[6]),
                        "exp": float(fields[9]),
                        "description": " ".join(fields[18:])
                    })
        df = pd.DataFrame(records)
        df.to_csv(out_csv, index=False)
        self.log_signal.emit(f"[✓] Parsed tblout saved to: {out_csv}")

    def parse_domtblout(self, domtbl_file, out_csv):
        self.log_signal.emit(f"[i] Parsing domtblout: {domtbl_file}")
        records = []
        with open(domtbl_file) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split()
                if len(fields) >= 23:
                    records.append({
                        "target_name": fields[0],
                        "target_accession": fields[1],
                        "query_name": fields[3],
                        "E-value": float(fields[6]),
                        "score": float(fields[7]),
                        "bias": float(fields[8]),
                        "domain_coord": f"{fields[17]}-{fields[18]}",
                        "description": " ".join(fields[22:])
                    })
        df = pd.DataFrame(records)
        df.to_csv(out_csv, index=False)
        self.log_signal.emit(f"[✓] Parsed domtblout saved to: {out_csv}")


# ====================== SCIENTIFIC INPUT WIDGET ======================
class ScientificInput(QWidget):
    def __init__(self, label, default_value, parent=None):
        super().__init__(parent)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        
        self.label = QLabel(label)
        layout.addWidget(self.label)
        
        self.coeff_input = QLineEdit()
        self.coeff_input.setFixedWidth(50)
        self.coeff_input.setText(str(default_value).split('e')[0])
        layout.addWidget(self.coeff_input)
        
        layout.addWidget(QLabel("× 10"))
        
        self.exponent_input = QLineEdit()
        self.exponent_input.setFixedWidth(50)
        
        # Try to extract exponent from scientific notation
        if 'e' in str(default_value):
            try:
                exp = str(default_value).split('e')[1]
                if exp.startswith('-'):
                    exp = exp[1:]
                self.exponent_input.setText(exp)
            except:
                self.exponent_input.setText("10")
        else:
            self.exponent_input.setText("0")
        
        layout.addWidget(self.exponent_input)
        
        layout.addWidget(QLabel("(negative exponent)"))
        
        # Set validators
        float_validator = QRegularExpressionValidator(QRegularExpression(r"[0-9]*\.?[0-9]+"))
        self.coeff_input.setValidator(float_validator)
        
        int_validator = QRegularExpressionValidator(QRegularExpression(r"[0-9]+"))
        self.exponent_input.setValidator(int_validator)
        
        # Add example
        self.example_label = QLabel("Example: 1.23e-10 → coeff: 1.23, exponent: 10")
        layout.addWidget(self.example_label)
        
        layout.addStretch()

    def value(self):
        try:
            coeff = float(self.coeff_input.text())
            exponent = int(self.exponent_input.text())
            return coeff * (10 ** -exponent)
        except:
            return None

    def set_value(self, value):
        if value == 0:
            self.coeff_input.setText("0")
            self.exponent_input.setText("0")
            return
        
        # Convert to scientific notation
        sci_notation = "{:.2e}".format(value)
        coeff, exp = sci_notation.split('e')
        exp = exp.replace('-', '')
        self.coeff_input.setText(coeff)
        self.exponent_input.setText(exp)


# ====================== MAIN GUI ======================
class HmmerGUI(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("HMMER Scanner")
        self.setGeometry(100, 100, 1000, 750)
        
        # Create menu bar
        self.create_menu_bar()
        
        # Central widget and main layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QVBoxLayout(central_widget)
        
        # Create tab widget
        tab_widget = QTabWidget()
        main_layout.addWidget(tab_widget)
        
        # Input tab
        input_tab = QWidget()
        tab_widget.addTab(input_tab, "Input")
        input_layout = QVBoxLayout(input_tab)
        
        # Input group
        input_group = QGroupBox("HMMER Configuration")
        input_layout.addWidget(input_group)
        form_layout = QVBoxLayout(input_group)
        
        # HMM Profile
        hmm_layout = QHBoxLayout()
        self.hmm_label = QLabel("HMM Profile:")
        self.hmm_input = QLineEdit()
        self.hmm_input.setPlaceholderText("/path/to/your/profiles.hmm")
        self.hmm_browse = QPushButton("Browse...")
        self.hmm_browse.clicked.connect(self.browse_hmm)
        hmm_layout.addWidget(self.hmm_label)
        hmm_layout.addWidget(self.hmm_input, 1)
        hmm_layout.addWidget(self.hmm_browse)
        form_layout.addLayout(hmm_layout)
        
        # FASTA File
        fasta_layout = QHBoxLayout()
        self.fasta_label = QLabel("FASTA File:")
        self.fasta_input = QLineEdit()
        self.fasta_input.setPlaceholderText("/path/to/your/query.fasta")
        self.fasta_browse = QPushButton("Browse...")
        self.fasta_browse.clicked.connect(self.browse_fasta)
        fasta_layout.addWidget(self.fasta_label)
        fasta_layout.addWidget(self.fasta_input, 1)
        fasta_layout.addWidget(self.fasta_browse)
        form_layout.addLayout(fasta_layout)
        
        # Output Directory
        output_dir_layout = QHBoxLayout()
        self.output_dir_label = QLabel("Output Directory:")
        self.output_dir_input = QLineEdit(os.getcwd())
        self.output_dir_browse = QPushButton("Browse...")
        self.output_dir_browse.clicked.connect(self.browse_output_dir)
        output_dir_layout.addWidget(self.output_dir_label)
        output_dir_layout.addWidget(self.output_dir_input, 1)
        output_dir_layout.addWidget(self.output_dir_browse)
        form_layout.addLayout(output_dir_layout)
        
        # Output Base
        output_layout = QHBoxLayout()
        self.output_label = QLabel("Output Base Name:")
        self.output_input = QLineEdit("hmmscan_results")
        output_layout.addWidget(self.output_label)
        output_layout.addWidget(self.output_input, 1)
        form_layout.addLayout(output_layout)
        
        # Thresholds group
        thresholds_group = QGroupBox("Threshold Settings")
        thresholds_layout = QVBoxLayout(thresholds_group)
        form_layout.addWidget(thresholds_group)
        
        # E-value
        self.evalue_widget = ScientificInput("Sequence E-value:", 1e-10)
        thresholds_layout.addWidget(self.evalue_widget)
        
        # Domain E-value
        self.dom_evalue_widget = ScientificInput("Domain E-value:", 1e-10)
        thresholds_layout.addWidget(self.dom_evalue_widget)
        
        # CPU
        cpu_layout = QHBoxLayout()
        self.cpu_label = QLabel("CPU Cores:")
        self.cpu_spin = QSpinBox()
        self.cpu_spin.setRange(1, os.cpu_count() or 4)
        self.cpu_spin.setValue(4)
        cpu_layout.addWidget(self.cpu_label)
        cpu_layout.addWidget(self.cpu_spin)
        cpu_layout.addStretch()
        thresholds_layout.addLayout(cpu_layout)
        
        # Extra Options
        extra_layout = QHBoxLayout()
        self.extra_label = QLabel("Extra Options:")
        self.extra_input = QLineEdit("--notextw")
        extra_layout.addWidget(self.extra_label)
        extra_layout.addWidget(self.extra_input, 1)
        form_layout.addLayout(extra_layout)
        
        # Run button
        self.run_button = QPushButton("Run HMMER Scan")
        self.run_button.setStyleSheet("""
            QPushButton {
                background-color: #4CAF50;
                color: white;
                font-weight: bold;
                padding: 10px;
                border-radius: 5px;
                font-size: 14px;
            }
            QPushButton:hover {
                background-color: #45a049;
            }
            QPushButton:disabled {
                background-color: #cccccc;
            }
        """)
        self.run_button.clicked.connect(self.run_scan)
        form_layout.addWidget(self.run_button)
        
        # Output tab
        output_tab = QWidget()
        tab_widget.addTab(output_tab, "Output")
        output_layout = QVBoxLayout(output_tab)
        
        # Output console
        self.output_console = QTextEdit()
        self.output_console.setReadOnly(True)
        self.output_console.setFont(QFont("Courier", 10))
        
        # Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        self.progress_bar.setFormat("Idle")
        self.progress_bar.setTextVisible(True)
        self.progress_bar.setStyleSheet("""
            QProgressBar {
                border: 1px solid #cccccc;
                border-radius: 3px;
                text-align: center;
            }
            QProgressBar::chunk {
                background-color: #4CAF50;
                width: 10px;
            }
        """)
        
        # Output files group
        files_group = QGroupBox("Output Files")
        files_layout = QVBoxLayout(files_group)
        
        self.tblout_label = QLabel("Sequence Hits (tblout):")
        self.tblout_label.setWordWrap(True)
        self.domtblout_label = QLabel("Domain Hits (domtblout):")
        self.domtblout_label.setWordWrap(True)
        self.output_dir_display = QLabel("Output Directory:")
        self.output_dir_display.setWordWrap(True)
        
        files_layout.addWidget(QLabel("Results will appear here after scan completes"))
        files_layout.addWidget(self.output_dir_display)
        files_layout.addWidget(self.tblout_label)
        files_layout.addWidget(self.domtblout_label)
        
        # Open folder button
        self.open_folder_button = QPushButton("Open Output Folder")
        self.open_folder_button.setStyleSheet("""
            QPushButton {
                background-color: #2196F3;
                color: white;
                padding: 5px;
                border-radius: 3px;
            }
            QPushButton:hover {
                background-color: #0b7dda;
            }
        """)
        self.open_folder_button.clicked.connect(self.open_output_folder)
        self.open_folder_button.setEnabled(False)
        files_layout.addWidget(self.open_folder_button)
        
        # Splitter for console and files
        splitter = QSplitter(Qt.Vertical)
        splitter.addWidget(self.output_console)
        splitter.addWidget(files_group)
        splitter.setSizes([500, 200])
        
        output_layout.addWidget(self.progress_bar)
        output_layout.addWidget(splitter, 1)
        
        # Status bar
        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)
        self.status_bar.showMessage("Ready to run HMMER scan")
        
        # Worker thread
        self.worker_thread = None
        self.worker = None
        
        # Set default theme
        self.set_dark_theme()

    def create_menu_bar(self):
        # Create menu bar
        menu_bar = QMenuBar()
        self.setMenuBar(menu_bar)
        
        # File menu
        file_menu = menu_bar.addMenu("&File")
        
        # New Project
        new_action = QAction("&New Project", self)
        new_action.setShortcut(QKeySequence.New)
        new_action.triggered.connect(self.new_project)
        file_menu.addAction(new_action)
        
        # Open Project
        open_action = QAction("&Open Project...", self)
        open_action.setShortcut(QKeySequence.Open)
        open_action.triggered.connect(self.open_project)
        file_menu.addAction(open_action)
        
        # Save Project
        save_action = QAction("&Save Project", self)
        save_action.setShortcut(QKeySequence.Save)
        save_action.triggered.connect(self.save_project)
        file_menu.addAction(save_action)
        
        file_menu.addSeparator()
        
        # Exit
        exit_action = QAction("E&xit", self)
        exit_action.setShortcut(QKeySequence.Quit)
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)
        
        # Edit menu
        edit_menu = menu_bar.addMenu("&Edit")
        
        # Copy
        copy_action = QAction("&Copy", self)
        copy_action.setShortcut(QKeySequence.Copy)
        copy_action.triggered.connect(self.copy_text)
        edit_menu.addAction(copy_action)
        
        # Clear
        clear_action = QAction("C&lear Console", self)
        clear_action.triggered.connect(self.clear_console)
        edit_menu.addAction(clear_action)
        
        # View menu
        view_menu = menu_bar.addMenu("&View")
        
        # Theme submenu
        theme_menu = view_menu.addMenu("&Theme")
        
        # Create theme action group
        self.theme_group = QActionGroup(self)
        self.theme_group.setExclusive(True)
        
        # Dark theme
        dark_theme_action = QAction("&Dark Theme", self)
        dark_theme_action.setCheckable(True)
        dark_theme_action.setChecked(True)
        dark_theme_action.triggered.connect(lambda: self.set_dark_theme())
        theme_menu.addAction(dark_theme_action)
        self.theme_group.addAction(dark_theme_action)
        
        # Light theme
        light_theme_action = QAction("&Light Theme", self)
        light_theme_action.setCheckable(True)
        light_theme_action.triggered.connect(lambda: self.set_light_theme())
        theme_menu.addAction(light_theme_action)
        self.theme_group.addAction(light_theme_action)
        
        # Help menu
        help_menu = menu_bar.addMenu("&Help")
        
        # Documentation
        docs_action = QAction("&Documentation", self)
        docs_action.triggered.connect(self.show_documentation)
        help_menu.addAction(docs_action)
        
        # About
        about_action = QAction("&About", self)
        about_action.triggered.connect(self.show_about)
        help_menu.addAction(about_action)

    def set_dark_theme(self):
        dark_palette = QPalette()
        dark_palette.setColor(QPalette.Window, QColor(53, 53, 53))
        dark_palette.setColor(QPalette.WindowText, Qt.white)
        dark_palette.setColor(QPalette.Base, QColor(25, 25, 25))
        dark_palette.setColor(QPalette.AlternateBase, QColor(53, 53, 53))
        dark_palette.setColor(QPalette.ToolTipBase, Qt.white)
        dark_palette.setColor(QPalette.ToolTipText, Qt.white)
        dark_palette.setColor(QPalette.Text, Qt.white)
        dark_palette.setColor(QPalette.Button, QColor(53, 53, 53))
        dark_palette.setColor(QPalette.ButtonText, Qt.white)
        dark_palette.setColor(QPalette.BrightText, Qt.red)
        dark_palette.setColor(QPalette.Link, QColor(42, 130, 218))
        dark_palette.setColor(QPalette.Highlight, QColor(42, 130, 218))
        dark_palette.setColor(QPalette.HighlightedText, Qt.black)
        
        QApplication.instance().setPalette(dark_palette)
        
        # Additional styling
        self.setStyleSheet("""
            QGroupBox {
                border: 1px solid #444;
                border-radius: 6px;
                margin-top: 10px;
                padding-top: 15px;
                font-weight: bold;
                background-color: #333;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                subcontrol-position: top left;
                padding: 0 5px;
                left: 10px;
                color: #eee;
            }
            QLabel {
                color: #ddd;
            }
            QLineEdit, QSpinBox, QDoubleSpinBox {
                background-color: #333;
                color: #fff;
                border: 1px solid #555;
                border-radius: 3px;
                padding: 3px;
            }
            QTextEdit {
                background-color: #252525;
                color: #ddd;
            }
            QTabWidget::pane {
                border: 1px solid #444;
                background: #333;
            }
            QTabBar::tab {
                background: #444;
                color: #ddd;
                padding: 8px;
                border-top-left-radius: 4px;
                border-top-right-radius: 4px;
            }
            QTabBar::tab:selected {
                background: #555;
                color: white;
            }
        """)
        self.status_bar.showMessage("Dark theme activated")

    def set_light_theme(self):
        light_palette = QApplication.style().standardPalette()
        QApplication.instance().setPalette(light_palette)
        
        # Additional styling
        self.setStyleSheet("""
            QGroupBox {
                border: 1px solid #ccc;
                border-radius: 6px;
                margin-top: 10px;
                padding-top: 15px;
                font-weight: bold;
                background-color: #f9f9f9;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                subcontrol-position: top left;
                padding: 0 5px;
                left: 10px;
                color: #333;
            }
            QLabel {
                color: #333;
            }
            QLineEdit, QSpinBox, QDoubleSpinBox {
                background-color: white;
                color: #000;
                border: 1px solid #ccc;
                border-radius: 3px;
                padding: 3px;
            }
            QTextEdit {
                background-color: white;
                color: #000;
            }
            QTabWidget::pane {
                border: 1px solid #ccc;
                background: #f0f0f0;
            }
            QTabBar::tab {
                background: #e0e0e0;
                color: #333;
                padding: 8px;
                border-top-left-radius: 4px;
                border-top-right-radius: 4px;
            }
            QTabBar::tab:selected {
                background: #f9f9f9;
                color: #000;
            }
        """)
        self.status_bar.showMessage("Light theme activated")

    def browse_hmm(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Select HMM Profile", "", "HMM Files (*.hmm)"
        )
        if file_path:
            self.hmm_input.setText(file_path)

    def browse_fasta(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Select FASTA File", "", "FASTA Files (*.fasta *.fa *.faa)"
        )
        if file_path:
            self.fasta_input.setText(file_path)

    def browse_output_dir(self):
        dir_path = QFileDialog.getExistingDirectory(
            self, "Select Output Directory", self.output_dir_input.text()
        )
        if dir_path:
            self.output_dir_input.setText(dir_path)

    def open_output_folder(self):
        output_dir = self.output_dir_input.text()
        if os.path.exists(output_dir):
            if sys.platform == "win32":
                os.startfile(output_dir)
            elif sys.platform == "darwin":
                subprocess.Popen(["open", output_dir])
            else:
                subprocess.Popen(["xdg-open", output_dir])
        else:
            QMessageBox.warning(self, "Directory Not Found", "Output directory does not exist")

    def run_scan(self):
        # Validate inputs
        if not self.hmm_input.text() or not os.path.exists(self.hmm_input.text()):
            QMessageBox.warning(self, "Missing Input", "Please select a valid HMM profile file")
            return
            
        if not self.fasta_input.text() or not os.path.exists(self.fasta_input.text()):
            QMessageBox.warning(self, "Missing Input", "Please select a valid FASTA file")
            return
            
        # Get E-values
        evalue = self.evalue_widget.value()
        dom_evalue = self.dom_evalue_widget.value()
        
        if evalue is None or dom_evalue is None:
            QMessageBox.warning(self, "Invalid Input", "Please enter valid E-values (e.g., coefficient: 1.0, exponent: 10)")
            return
            
        # Get output directory
        output_dir = self.output_dir_input.text()
        if not output_dir:
            output_dir = os.getcwd()
            
        # Create output directory if it doesn't exist
        try:
            os.makedirs(output_dir, exist_ok=True)
        except Exception as e:
            QMessageBox.critical(self, "Directory Error", f"Could not create output directory:\n\n{str(e)}")
            return
            
        # Disable UI during processing
        self.run_button.setEnabled(False)
        self.run_button.setText("Processing...")
        self.progress_bar.setValue(0)
        self.progress_bar.setFormat("Starting...")
        self.output_console.clear()
        
        # Prepare config
        config = {
            'hmm_profile': self.hmm_input.text(),
            'fasta_file': self.fasta_input.text(),
            'output_dir': output_dir,
            'output_base': self.output_input.text(),
            'evalue': evalue,
            'dom_evalue': dom_evalue,
            'cpu': self.cpu_spin.value(),
            'extra_options': self.extra_input.text().split()
        }
        
        # Create worker and thread
        self.worker = HmmerWorker(config)
        self.worker_thread = QThread()
        self.worker.moveToThread(self.worker_thread)
        
        # Connect signals
        self.worker.log_signal.connect(self.log_message)
        self.worker.progress_signal.connect(self.update_progress)
        self.worker.finished_signal.connect(self.scan_finished)
        self.worker.error_signal.connect(self.scan_error)
        self.worker.file_output_signal.connect(self.update_output_files)
        
        # Start thread
        self.worker_thread.started.connect(self.worker.run)
        self.worker_thread.start()
        
        self.status_bar.showMessage("HMMER scan in progress...")

    def log_message(self, message):
        self.output_console.append(message)
        self.output_console.moveCursor(QTextCursor.End)

    def update_progress(self, value):
        self.progress_bar.setValue(value)
        self.progress_bar.setFormat(f"Processing: {value}%")

    def update_output_files(self, file_type, path):
        if file_type == "tblout":
            self.tblout_label.setText(f"Sequence Hits (tblout): {path}")
        elif file_type == "domtblout":
            self.domtblout_label.setText(f"Domain Hits (domtblout): {path}")
        elif file_type == "stdout":
            pass  # We already log this
        elif file_type == "output_dir":
            self.output_dir_display.setText(f"Output Directory: {path}")
            self.open_folder_button.setEnabled(True)

    def scan_finished(self):
        self.worker_thread.quit()
        self.worker_thread.wait()
        
        self.run_button.setEnabled(True)
        self.run_button.setText("Run HMMER Scan")
        self.progress_bar.setFormat("Completed successfully")
        self.status_bar.showMessage("HMMER scan completed successfully")
        
        # Update output directory display
        output_dir = self.output_dir_input.text()
        self.output_dir_display.setText(f"Output Directory: {output_dir}")
        self.open_folder_button.setEnabled(True)
        
        QMessageBox.information(self, "Success", "HMMER scan completed successfully!")

    def scan_error(self, message):
        self.worker_thread.quit()
        self.worker_thread.wait()
        
        self.run_button.setEnabled(True)
        self.run_button.setText("Run HMMER Scan")
        self.progress_bar.setFormat("Error occurred")
        
        # Highlight error in red
        self.log_message(f"\n[!] ERROR: {message}")
        
        self.status_bar.showMessage(f"Error: {message}")
        QMessageBox.critical(self, "Error", f"HMMER scan failed:\n\n{message}")
    
    # Menu actions
    def new_project(self):
        # Reset all fields to default
        self.hmm_input.clear()
        self.fasta_input.clear()
        self.output_dir_input.setText(os.getcwd())
        self.output_input.setText("hmmscan_results")
        self.evalue_widget.set_value(1e-10)
        self.dom_evalue_widget.set_value(1e-10)
        self.cpu_spin.setValue(4)
        self.extra_input.setText("--notextw")
        self.output_console.clear()
        self.status_bar.showMessage("New project created")
        
    def open_project(self):
        # Placeholder for actual implementation
        QMessageBox.information(self, "Open Project", "Project opening functionality would be implemented here")
        
    def save_project(self):
        # Placeholder for actual implementation
        QMessageBox.information(self, "Save Project", "Project saving functionality would be implemented here")
        
    def copy_text(self):
        if self.output_console.hasFocus():
            self.output_console.copy()
        else:
            # Copy from other focused widgets
            focused = QApplication.focusWidget()
            if isinstance(focused, QLineEdit):
                focused.copy()
                
    def clear_console(self):
        self.output_console.clear()
        self.status_bar.showMessage("Console cleared")
        
    def show_documentation(self):
        QMessageBox.information(self, "Documentation", 
            "HMMER Scanner Documentation\n\n"
            "1. Select HMM profile and FASTA file\n"
            "2. Configure output settings and thresholds\n"
            "3. Click 'Run HMMER Scan' to start analysis\n"
            "4. View results in the Output tab\n\n"
            "For more information, visit the HMMER documentation: "
            "http://hmmer.org/documentation.html")
        
    def show_about(self):
        QMessageBox.about(self, "About HMMER Scanner",
            "HMMER Scanner\nVersion 1.0\n\n"
            "A graphical interface for running HMMER protein sequence analysis.\n\n"
            "Developed using PySide6 and HMMER tools.\n\n"
            "© 2023 BioInformatics Tools")


if __name__ == "__main__":
    app = QApplication(sys.argv)
    
    # Set application style
    app.setStyle("Fusion")
    
    window = HmmerGUI()
    window.show()
    sys.exit(app.exec())