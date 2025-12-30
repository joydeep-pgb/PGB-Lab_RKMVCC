import sys
import os
from PySide6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                             QHBoxLayout, QLabel, QLineEdit, QPushButton, 
                             QFileDialog, QSpinBox, QTextEdit, QGroupBox, 
                             QTabWidget, QProgressBar, QComboBox, QDoubleSpinBox,
                             QStyleFactory, QMessageBox)
from PySide6.QtCore import QProcess, Qt, QSettings
from PySide6.QtGui import QFont, QTextCursor

class EggNogPro(QMainWindow):
    def __init__(self):
        super().__init__()
        
        # Application Config
        self.setWindowTitle("eggNOG-mapper Professional")
        self.resize(900, 700)
        
        # Persistent Settings (Saves inputs between sessions)
        self.settings = QSettings("BioTools", "EggNogMapperPro")
        
        # State
        self.process = None
        self.is_running = False

        # Initialization
        self.setup_style()
        self.init_ui()
        self.load_settings()

    def setup_style(self):
        """Applies a professional dark theme."""
        QApplication.setStyle(QStyleFactory.create("Fusion"))
        self.setStyleSheet("""
            QMainWindow { background-color: #2b2b2b; color: #ffffff; }
            QGroupBox { 
                font-weight: bold; border: 1px solid #555; margin-top: 10px; padding-top: 15px; border-radius: 5px; color: #ddd;
            }
            QGroupBox::title { subcontrol-origin: margin; left: 10px; padding: 0 5px; }
            QLabel { color: #e0e0e0; }
            QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox { 
                background-color: #3b3b3b; color: #fff; border: 1px solid #555; padding: 5px; border-radius: 3px;
            }
            QPushButton {
                background-color: #3c3c3c; color: white; border: 1px solid #555; padding: 6px 12px; border-radius: 4px;
            }
            QPushButton:hover { background-color: #4c4c4c; }
            QPushButton:pressed { background-color: #2c2c2c; }
            QTabWidget::pane { border: 1px solid #444; }
            QTabBar::tab { background: #2b2b2b; color: #aaa; padding: 8px 20px; border: 1px solid #444; border-bottom: none; }
            QTabBar::tab:selected { background: #3b3b3b; color: #fff; border-bottom: 2px solid #2b78e4; }
            QProgressBar { text-align: center; border: 1px solid #555; border-radius: 3px; }
            QProgressBar::chunk { background-color: #2b78e4; }
        """)

    def init_ui(self):
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QVBoxLayout(central_widget)

        # --- Tabs ---
        self.tabs = QTabWidget()
        main_layout.addWidget(self.tabs)

        # Tab 1: Run (Basic Input/Output)
        self.tab_run = QWidget()
        self.setup_run_tab()
        self.tabs.addTab(self.tab_run, "Run Configuration")

        # Tab 2: Advanced Settings
        self.tab_advanced = QWidget()
        self.setup_advanced_tab()
        self.tabs.addTab(self.tab_advanced, "Advanced Options")

        # --- Logging Area ---
        log_group = QGroupBox("Execution Log")
        log_layout = QVBoxLayout(log_group)
        self.log_window = QTextEdit()
        self.log_window.setReadOnly(True)
        self.log_window.setFont(QFont("Consolas", 10))
        self.log_window.setStyleSheet("background-color: #1e1e1e; color: #d4d4d4; border: none;")
        log_layout.addWidget(self.log_window)
        main_layout.addWidget(log_group, stretch=1)

        # --- Control Bar ---
        control_layout = QHBoxLayout()
        
        self.status_label = QLabel("Ready")
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 0) # Indeterminate mode default
        self.progress_bar.setTextVisible(False)
        self.progress_bar.setFixedSize(150, 15)
        self.progress_bar.hide()

        self.btn_run = QPushButton("Start Annotation")
        self.btn_run.setFixedHeight(40)
        self.btn_run.setStyleSheet("background-color: #2b78e4; font-weight: bold;")
        self.btn_run.clicked.connect(self.toggle_process)

        self.btn_stop = QPushButton("Stop")
        self.btn_stop.setFixedHeight(40)
        self.btn_stop.setStyleSheet("background-color: #c42b2b; font-weight: bold;")
        self.btn_stop.setEnabled(False)
        self.btn_stop.clicked.connect(self.kill_process)

        control_layout.addWidget(self.status_label)
        control_layout.addWidget(self.progress_bar)
        control_layout.addStretch()
        control_layout.addWidget(self.btn_stop)
        control_layout.addWidget(self.btn_run)
        
        main_layout.addLayout(control_layout)

        # Process Setup
        self.process = QProcess(self)
        self.process.readyReadStandardOutput.connect(self.handle_stdout)
        self.process.readyReadStandardError.connect(self.handle_stderr)
        self.process.finished.connect(self.process_finished)

    def setup_run_tab(self):
        layout = QVBoxLayout(self.tab_run)
        
        # Files Group
        file_group = QGroupBox("Required Inputs")
        f_layout = QVBoxLayout(file_group)
        
        self.input_file = self.create_browser_row(f_layout, "Input FASTA:", "Select Input", True)
        self.db_dir = self.create_browser_row(f_layout, "Database Path:", "Select DB Dir", False)
        self.output_prefix = self.create_text_row(f_layout, "Output Prefix:")
        
        layout.addWidget(file_group)
        layout.addStretch()

    def setup_advanced_tab(self):
        layout = QVBoxLayout(self.tab_advanced)
        
        # Performance Group
        perf_group = QGroupBox("Performance")
        p_layout = QHBoxLayout(perf_group)
        
        p_layout.addWidget(QLabel("CPU Threads:"))
        self.cpu_spin = QSpinBox()
        self.cpu_spin.setRange(1, 128)
        self.cpu_spin.setValue(8)
        self.cpu_spin.setToolTip("Number of CPU threads to use.")
        p_layout.addWidget(self.cpu_spin)
        
        p_layout.addWidget(QLabel("   Mode:"))
        self.mode_combo = QComboBox()
        self.mode_combo.addItems(["diamond", "mmseqs"])
        self.mode_combo.setToolTip("Search tool to use.")
        p_layout.addWidget(self.mode_combo)
        p_layout.addStretch()
        
        # Filtering Group
        filter_group = QGroupBox("Filtering Parameters")
        fl_layout = QHBoxLayout(filter_group)
        
        # Updated Defaults per Manual: 0.001 E-value, 0 Score
        fl_layout.addWidget(QLabel("Min E-value:"))
        self.evalue_spin = QDoubleSpinBox()
        self.evalue_spin.setDecimals(5)
        self.evalue_spin.setRange(0, 10)
        self.evalue_spin.setValue(0.001)
        self.evalue_spin.setToolTip("Default: 0.001")
        fl_layout.addWidget(self.evalue_spin)
        
        fl_layout.addWidget(QLabel("   Min Bit Score:"))
        self.score_spin = QDoubleSpinBox()
        self.score_spin.setRange(0, 99999) # Allow high scores
        self.score_spin.setValue(0)
        self.score_spin.setToolTip("Default: 0 (No threshold)")
        fl_layout.addWidget(self.score_spin)
        
        fl_layout.addStretch()

        layout.addWidget(perf_group)
        layout.addWidget(filter_group)
        layout.addStretch()

    # --- Helper Widgets ---
    def create_browser_row(self, layout, label_text, btn_text, is_file):
        row = QHBoxLayout()
        lbl = QLabel(label_text)
        lbl.setFixedWidth(100)
        row.addWidget(lbl)
        
        line_edit = QLineEdit()
        row.addWidget(line_edit)
        
        btn = QPushButton(btn_text)
        if is_file:
            btn.clicked.connect(lambda: self.browse_file(line_edit))
        else:
            btn.clicked.connect(lambda: self.browse_dir(line_edit))
        row.addWidget(btn)
        
        layout.addLayout(row)
        return line_edit

    def create_text_row(self, layout, label_text):
        row = QHBoxLayout()
        lbl = QLabel(label_text)
        lbl.setFixedWidth(100)
        row.addWidget(lbl)
        line_edit = QLineEdit()
        row.addWidget(line_edit)
        layout.addLayout(row)
        return line_edit

    # --- Logic ---
    def browse_file(self, line_edit):
        f, _ = QFileDialog.getOpenFileName(self, "Select FASTA", "", "FASTA (*.fa *.faa *.fasta);;All (*.*)")
        if f: line_edit.setText(f)

    def browse_dir(self, line_edit):
        d = QFileDialog.getExistingDirectory(self, "Select Database Directory")
        if d: line_edit.setText(d)

    def toggle_process(self):
        if self.is_running:
            return 
        
        # Validation
        if not self.input_file.text() or not self.db_dir.text():
            QMessageBox.warning(self, "Missing Input", "Please provide both Input File and Database Directory.")
            return

        # Prepare UI
        self.is_running = True
        self.btn_run.setEnabled(False)
        self.btn_stop.setEnabled(True)
        self.progress_bar.show()
        self.status_label.setText("Running...")
        self.log_window.clear()
        self.save_settings()

        # Build Command
        program = "emapper.py" 
        
        # Construct arguments
        args = [
            "-i", self.input_file.text(),
            "--data_dir", self.db_dir.text(),
            "-o", self.output_prefix.text() or "results",
            "-m", self.mode_combo.currentText(),
            "--cpu", str(self.cpu_spin.value()),
            "--evalue", str(self.evalue_spin.value()),
            "--score", str(self.score_spin.value()),
            "--override"
        ]

        self.log_window.append(f"<span style='color: #2b78e4;'><b>Command:</b> {program} {' '.join(args)}</span><br>")
        
        # Start Process
        self.process.start(program, args)

    def kill_process(self):
        if self.process.state() == QProcess.Running:
            self.process.kill()
            self.log_window.append("<br><span style='color: red;'><b>Process Terminated by User</b></span>")

    def handle_stdout(self):
        data = self.process.readAllStandardOutput().data().decode()
        self.log_window.moveCursor(QTextCursor.End)
        self.log_window.insertPlainText(data)
        self.log_window.moveCursor(QTextCursor.End)

    def handle_stderr(self):
        data = self.process.readAllStandardError().data().decode()
        self.log_window.moveCursor(QTextCursor.End)
        # Highlight stderr in orange to differentiate from normal logs
        self.log_window.insertHtml(f"<span style='color: #ff9d00;'>{data}</span>") 
        self.log_window.moveCursor(QTextCursor.End)

    def process_finished(self):
        self.is_running = False
        self.btn_run.setEnabled(True)
        self.btn_stop.setEnabled(False)
        self.progress_bar.hide()
        
        exit_code = self.process.exitCode()
        if exit_code == 0:
            self.status_label.setText("Finished Successfully")
            self.log_window.append("<br><span style='color: #00ff00;'><b>--- Task Completed ---</b></span>")
        else:
            self.status_label.setText(f"Failed (Code {exit_code})")

    # --- Persistence ---
    def load_settings(self):
        self.input_file.setText(self.settings.value("input_file", ""))
        self.db_dir.setText(self.settings.value("db_dir", ""))
        self.output_prefix.setText(self.settings.value("out_prefix", "results"))
        self.cpu_spin.setValue(int(self.settings.value("cpu", 8)))

    def save_settings(self):
        self.settings.setValue("input_file", self.input_file.text())
        self.settings.setValue("db_dir", self.db_dir.text())
        self.settings.setValue("out_prefix", self.output_prefix.text())
        self.settings.setValue("cpu", self.cpu_spin.value())

if __name__ == "__main__":
    app = QApplication(sys.argv)
    gui = EggNogPro()
    gui.show()
    sys.exit(app.exec())